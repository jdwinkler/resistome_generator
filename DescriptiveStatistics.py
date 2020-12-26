from collections import defaultdict
from resistome.graphics import visualization as glib
from resistome import constants
from resistome.utils import database_utils, network_utils
import os
import numpy
import scipy
import scipy.stats
import scipy.misc
from resistome.sql import sql_interface


def interaction_matrix(source_data,
                       rows,
                       genome_size,
                       sig_threshold=0.0001,
                       ignore_tuples=set(),
                       remove_empty_rows=False):
    """
    
    Computes a fisher's exact test to evaluate enrichment of the intersection between source_data[x] ^ source_data[y].
    
    Intended to be used for interaction detection between different phenotypes; possible data for source_data[x]:
    
    *gene names
    *metabolite names
    *go processes (biological function, molecular role, structural component)
    
    This function does not presume any particular data input and should work without modification for other future
    data additions.
    
    :param source_data: map of x := list of data for each x key
    :param rows: iterable (should be a list), keys in source_data
    :param genome_size: integer, number of elements in "genome" (sometimes literally, sometimes number of metabolites
                        or GO numbers; user responsible for providing correct data
    :param sig_threshold: float/double, minimum P-value by fisher's exact test to consider significant
    :param ignore_tuples: set of tuples, interaction (x, y) tuples to ignore (x, y are strings)
    :param remove_empty_rows: boolean, include rows that have no significant interactions
    :return: 
    
    """

    a_mat = defaultdict(dict)

    significant_tuples = []

    for cond_x in rows:

        x_set = set(source_data[cond_x])
        x_len = len(x_set)

        for cond_y in rows:

            if cond_x == cond_y:
                a_mat[cond_x][cond_y] = (0.0, 0.0)
                continue

            if (cond_x, cond_y) in ignore_tuples:
                a_mat[cond_x][cond_y] = (0.0, 0.0)
                a_mat[cond_y][cond_x] = (0.0, 0.0)
                continue

            if cond_x in source_data and cond_y in source_data:

                if cond_y in a_mat[cond_x]:
                    continue

                y_set = set(source_data[cond_y])
                y_len = len(set(source_data[cond_y]))
                overlap = len(x_set & y_set)
                expected_number = float(x_len) * float(y_len) / float(genome_size)

                if expected_number == 0:
                    a_mat[cond_x][cond_y] = (0.0, 0.0)
                    a_mat[cond_y][cond_x] = (0.0, 0.0)
                    continue

                contingency_table = [[overlap, len(x_set - y_set)],
                                     [len(y_set - x_set), genome_size - overlap - len(x_set - y_set) - len(y_set - x_set)]]

                all_above_min = True
                for row in contingency_table:
                    for col in row:
                        if col < 10:
                            all_above_min = False
                            break

                if all_above_min:
                    (oddsratio, probability, _, _,) = scipy.stats.chi2_contingency(contingency_table)
                else:
                    (oddsratio, probability) = scipy.stats.fisher_exact(contingency_table)

                if numpy.isnan(oddsratio) or numpy.isinf(oddsratio):
                    oddsratio = 50

                # basically, these values are very, very close to P = 0,
                # but still have an exponent of something like -30
                # if we represent these values like so, they remain very
                # significant but without the burden of make an arbitrary
                # choice about how to represent them.
                # still weird though.
                if probability == 0:
                    a_mat[cond_x][cond_y] = (0.0, 0.0)
                    a_mat[cond_y][cond_x] = (0.0, 0.0)
                else:
                    a_mat[cond_x][cond_y] = (abs(numpy.log10(probability)), oddsratio)
                    a_mat[cond_y][cond_x] = (abs(numpy.log10(probability)), oddsratio)

                if probability < sig_threshold:
                    significant_tuples.append((cond_x, cond_y, probability))

    # now pruned.
    rows = sorted(a_mat.keys(), reverse=True)

    if remove_empty_rows:
        keys_to_run = set()
        for x in rows:
            for y in a_mat[x]:
                if a_mat[x][y][0] > abs(numpy.log10(sig_threshold)):
                    keys_to_run.add(x)
                    break
        rows = keys_to_run

    data = []
    for cond_x in rows:

        columns = []
        for cond_y in rows:
            if cond_y in a_mat[cond_x]:
                (pvalue, odds) = a_mat[cond_x][cond_y]
                if pvalue > abs(numpy.log10(sig_threshold)):
                    columns.append(pvalue)
                else:
                    columns.append(0)
            else:
                columns.append(0)

        data.append(columns)

    matrix = numpy.matrix(data)
    return matrix, rows, significant_tuples


def association_btw_phenotypes(genome_size,
                               identifier,
                               types_of_interest=None,
                               output_prefix='',
                               ignore_threshold=0,
                               use_tags=True):

    def filter_correlated_tuples(tagset, paper_doi_dict, ignore_threshold):

        r_ignore_set = set()

        for tag_x in tagset:
            for tag_y in tagset:

                if tag_x == tag_y:
                    continue
                else:
                    x_dois = paper_doi_dict[tag_x]
                    y_dois = paper_doi_dict[tag_y]

                    if len(x_dois ^ y_dois) < ignore_threshold:
                        r_ignore_set.add((tag_x, tag_y))
                        r_ignore_set.add((tag_y, tag_x))

        return r_ignore_set

    if identifier not in {'gene', 'go', 'metabolite'}:
        raise AssertionError('Unknown id requested: %s' % identifier)

    sensitive_tag_overlap = defaultdict(set)
    resistant_tag_overlap = defaultdict(set)
    r_paper_counter = defaultdict(set)
    s_paper_counter = defaultdict(set)

    genotypes = sql_interface.get_mutant_genotypes_as_features(identifier)

    for result in genotypes:

        if use_tags:
            sql_key = 'pclasses'
        else:
            sql_key = 'phenotypes'

        if types_of_interest is not None and set(result['pclasses']).isdisjoint(set(types_of_interest)):
            continue

        for phenotype, ptype in zip(result[sql_key], result['types']):

            phenotype = phenotype.upper()
            features = result[identifier]

            if identifier == 'go':
                features = set(features)

            if ptype == 'S':
                sensitive_tag_overlap[phenotype].update(features)
                s_paper_counter[phenotype].add(result['paper_id'])
            else:
                resistant_tag_overlap[phenotype].update(features)
                r_paper_counter[phenotype].add(result['paper_id'])

    r_tuples = filter_correlated_tuples(resistant_tag_overlap.keys(),
                                        r_paper_counter, ignore_threshold)
    s_tuples = filter_correlated_tuples(sensitive_tag_overlap.keys(),
                                        s_paper_counter, ignore_threshold)

    rows = list(set(resistant_tag_overlap.keys()) | set(sensitive_tag_overlap.keys()))
    converted_ids = sql_interface.convert_full_identifiers_to_abbreviations('phenotype', rows)
    id_pairs = sorted([(phenotype,y) for phenotype,y in zip(rows, converted_ids)], key=lambda x: x[1], reverse=True)
    rows = [phenotype[0] for phenotype in id_pairs]
    display_rows = [phenotype[1] for phenotype in id_pairs]

    matrix, rows, resistant_sig_tuple = interaction_matrix(
        resistant_tag_overlap, rows,
        genome_size=genome_size,
        sig_threshold=1e-5,
        ignore_tuples=r_tuples)

    glib.generate_heatmap(numpy.array(matrix),
                          display_rows,
                          'log10[P-value]',
                          os.path.join(constants.OUTPUT_DIR,
                                       output_prefix + 'Resistance Interaction HeatMap.pdf'),
                          cmap='Greens')

    matrix, rows, sensitive_sig_tuple = interaction_matrix(
        sensitive_tag_overlap, rows,
        genome_size=genome_size,
        ignore_tuples=s_tuples)

    glib.generate_heatmap(numpy.array(matrix),
                          display_rows,
                          'log10[P-value]',
                          os.path.join(constants.OUTPUT_DIR,
                                       output_prefix + 'Sensitive Interaction HeatMap.pdf'),
                          cmap='Greens')

    return resistant_sig_tuple, sensitive_sig_tuple, resistant_tag_overlap, sensitive_tag_overlap


def map_genes_by_tag(genes_to_include, filename, tag_threshold=0, aggregate_threshold=5, binary=True):

    gene_tag_dict = defaultdict(dict)
    tag_list = set()

    seen_before = set()

    distribution_data = defaultdict(int)

    results = sql_interface.get_mutant_genotypes()
    for r in results:
        phenotype_classes = r['pclasses']
        phenotype_types = r['types']
        std_mg1655_names = r['mg1655_names']
        paper_id = r['paper_id']

        for std_mg1655_name in std_mg1655_names:
            for (phenotype_class, phenotype_type) in zip(phenotype_classes, phenotype_types):
                if (std_mg1655_name, phenotype_class, phenotype_type, paper_id) in seen_before:
                    continue
                else:
                    seen_before.add((std_mg1655_name, phenotype_class, phenotype_type, paper_id))

                if phenotype_class not in gene_tag_dict[std_mg1655_name]:
                    gene_tag_dict[std_mg1655_name][phenotype_class] = 1
                else:
                    gene_tag_dict[std_mg1655_name][phenotype_class] += 1

                tag_list.add(phenotype_class)

    genes_to_process = set()

    if len(genes_to_include) == 0:
        genes_to_include = gene_tag_dict.keys()

    for gene in genes_to_include:
        psum = 0
        distribution_data[gene] = sum(gene_tag_dict[gene].values())
        for phenotype_class in gene_tag_dict[gene]:
            psum += gene_tag_dict[gene][phenotype_class]
            if gene_tag_dict[gene][phenotype_class] > tag_threshold:
                genes_to_process.add(gene)
        if psum < aggregate_threshold and gene in genes_to_process:
            genes_to_process.remove(gene)

    gene_list = sorted(genes_to_process, reverse=True)
    tag_list = sorted(list(tag_list))

    output_tags = sql_interface.convert_full_identifiers_to_abbreviations('phenotype',
                                                                          tag_list)
    sorted_tags = sorted([(x, y) for x, y in zip(tag_list, output_tags)],
                         key=lambda x: x[1])

    tag_list = [x[0] for x in sorted_tags]
    output_tags = [x[1] for x in sorted_tags]

    matrix = []
    output_gene_names = []

    common_gene_dict = sql_interface.get_gene_common_names(genes_to_process)

    for gene in gene_list:

        if gene in common_gene_dict:
            output_gene_names.append(common_gene_dict[gene] + ' (%s)' % gene)
        else:
            output_gene_names.append(gene)

        vector = []
        for tag in tag_list:
            if tag in gene_tag_dict[gene] and gene_tag_dict[gene][tag] > tag_threshold:
                vector.append(1 if binary else min(gene_tag_dict[gene][tag], 5))
            else:
                vector.append(0)
        matrix.append(vector)

    output_gene_names = [x.upper() for x in output_gene_names]

    glib.generate_heatmap_asymmetric(numpy.array(matrix),
                                     output_tags,
                                     output_gene_names,
                                     'Genotype-Phenotype Occurrences',
                                     filename,
                                     cmap='Greys' if binary else 'Reds',
                                     plot_color_bar=True,
                                     dim_override=(8, max(0.279 * len(gene_list), 12)),
                                     x_label_rotation='vertical',
                                     y_label_rotation='horizontal')

    return distribution_data


def vital_statistics():

    def dict_to_paired_array(dictionary):

        tphen = []

        for key in dictionary:
            tphen.append((key, dictionary[key]))

        tphen = sorted(tphen, key=lambda x: x[1], reverse=True)

        tolerance_keys = [item[0] for item in tphen]
        tolerance_vals = [item[1] for item in tphen]

        return tolerance_keys, tolerance_vals

    year_papers, \
    year_models, \
    year_mutation, \
    paper_count, \
    specific_phenotypes, \
    method_counts, \
    num_ecoli_designs, \
    num_ecoli_genes, \
    design_results, \
    gene_usage, \
    tolerance_phenotypes, \
    mutantCounts, \
    mutationCounts, \
    mutation_counts, = sql_interface.get_table_statistics()


    design_values = []
    for r in design_results:
        design_values.append(r['count'])
    design_values = list(filter(lambda x: x > 1, design_values))
    more_than_fraction = len(design_values)

    gene_output = []

    for (gene, count) in gene_usage.items():
        gene_output.append((gene, count))

    go = sorted(gene_output, key=lambda x: x[1], reverse=True)
    go_values = [x[1] for x in go]

    topX = 25

    fraction = float(sum(go_values[0:topX - 1])) / float(sum(go_values))

    go = go[0:topX - 1]

    glib.bargraph([x[0] for x in go],
                  [x[1] for x in go],
                         'Mutated Genes',
                         'Mutation Counts',
                  os.path.join(constants.OUTPUT_DIR,
                                      'Mutated Gene Counts.pdf'),
                  rotation='vertical')

    glib.histogram(gene_usage.values(),
                   'Distribution of Gene Usage Counts',
                   'Counts',
                   os.path.join(constants.OUTPUT_DIR,
                                'Gene Usage Distribution.pdf'),
                   use_median_filter=True)

    tolerance_keys, tolerance_vals = dict_to_paired_array(tolerance_phenotypes)
    tolerance_keys = sql_interface.convert_full_identifiers_to_abbreviations('phenotype', tolerance_keys)
    glib.bargraph(tolerance_keys,
                  tolerance_vals,
                         'Tolerance Phenotype',
                         'Mutant Counts',
                  os.path.join(constants.OUTPUT_DIR,
                                      'Tolerance Counts.pdf'),
                  rotation='vertical')

    method_keys, method_vals = dict_to_paired_array(method_counts)
    glib.bargraph(method_keys,
                  method_vals,
                         'Engineering Method',
                         'Mutant Counts',
                  os.path.join(constants.OUTPUT_DIR,
                                      'Method Counts.pdf'),
                  rotation='vertical',
                  mapper=None)

    mutation_keys, mutation_vals = dict_to_paired_array(mutation_counts)
    mutation_keys = sql_interface.convert_full_identifiers_to_abbreviations('mutations', mutation_keys)
    glib.bargraph(mutation_keys,
                  mutation_vals,
                         'Mutation Type',
                         'Counts',
                  os.path.join(constants.OUTPUT_DIR,
                                      'Mutation Counts.pdf'),
                  rotation='vertical')

    paper_by_year = []
    mutations_by_year = []

    chronological = sorted(year_papers.keys())

    for year in chronological:
        paper_by_year.append((year, year_papers[year]))
        mutations_by_year.append((year, numpy.max(year_models[year])))

    # papers/mutations per year
    paper_values = [item[1] for item in paper_by_year]
    mutation_values = [item[1] for item in mutations_by_year]
    glib.linebar(chronological,
                 paper_values,
                 mutation_values,
                 'Publication Year',
                 'Resistome Papers',
                 'Mutants/Paper',
                 os.path.join(constants.OUTPUT_DIR, 'Papers and Mutants Per Year.pdf'))

    mutants_per_paper_avg = numpy.mean(mutantCounts)
    mutants_per_paper_std = numpy.std(mutantCounts)

    mutations_per_design_avg = numpy.mean(mutationCounts)
    mutations_per_design_std = numpy.std(mutationCounts)

    vital_stats = []

    vital_stats.append('The top %i genes represent %f fraction of all mutations (ratio of actual to expected is %f)' % \
          (topX, fraction,
           fraction / (float(topX) / float(len(gene_usage.keys())))))

    # vital_stats.append('Number of total mutants per paper: %f' % (total_per_paper_med))
    vital_stats.append('Number of recorded mutants per paper: %f (%f)' % (mutants_per_paper_avg, mutants_per_paper_std))
    vital_stats.append('Number of mutations per design %f (%f)' % (mutations_per_design_avg, mutations_per_design_std))
    vital_stats.append('Number of unique papers %i' % paper_count)
    vital_stats.append('Number of E. coli designs %i' % num_ecoli_designs)
    vital_stats.append('Number of unique phenotypes (both resistance + sensitivity) %i' % len(specific_phenotypes))
    vital_stats.append('Number of phenotype categories: %i' % len(tolerance_phenotypes.keys()))
    vital_stats.append('Number of E. coli genes: %i' % num_ecoli_genes)
    vital_stats.append('Number of unique genes %i' % len(gene_usage.keys()))
    vital_stats.append('Number of designs with more than 1 mutation: %i' % more_than_fraction)

    fhandle = open(os.path.join(constants.OUTPUT_DIR, 'Vital Statistics Report.txt'), 'w')

    for line in vital_stats:
        fhandle.write(line + '\n')

    fhandle.close()


def basic_information_figures():

    jkeys, jvalues, meth_keys, meth_values, mut_keys, mut_values = sql_interface.journal_method_design_distributions()

    glib.bargraph(jkeys, jvalues, 'Journals', 'Count',
                  os.path.join(constants.OUTPUT_DIR, 'Journal Distribution.pdf'),
                  rotation='vertical')

    # methods used for mutant design figure
    glib.bargraph(meth_keys, meth_values, 'Design Methods', 'Count',
                  os.path.join(constants.OUTPUT_DIR, 'Method Distribution.pdf'),
                  rotation='vertical')

    # mutations made figure
    glib.bargraph(mut_keys, mut_values, 'Mutation Types', 'Count',
                  os.path.join(constants.OUTPUT_DIR, 'Mutation Distribution.pdf'),
                  rotation='vertical')


def tag_distribution(output_file_name, top=9):

    frequency, category_counter = sql_interface.get_tag_statistics()

    phenotypes = sql_interface.get_phenotype_classes()
    converted_phenotypes = sql_interface.convert_full_identifiers_to_abbreviations('phenotype', keys=phenotypes)
    name_dict = dict()

    for p, con_p in zip(phenotypes, converted_phenotypes):
        name_dict[p.lower()] = con_p

    glib.generate_stacked(frequency,
                          category_counter,
                          'Year',
                          'Proportion',
                          os.path.join(constants.OUTPUT_DIR,
                                       output_file_name),
                          top=top,
                          cmap='tab20',
                          name_mapper=name_dict)


def random_versus_designed_comparison():

    paper_doe_types = set()
    paper_doe_types.update(sql_interface.get_paper_doe_types())

    non_random_types = paper_doe_types - {'random'}
    random_types = {'random'}

    random_types_values, random_gene_counts = sql_interface.get_mutation_counts_by_doe_type(types_to_include=random_types)

    nonrandom_types_values, nr_gene_counts = sql_interface.get_mutation_counts_by_doe_type(types_to_include=non_random_types)

    text_output = []

    text_output.append('Average number of mutational types per random, non-random studies: %f, %f, P-value (t-test, two-tailed) = %g' \
          % (numpy.average(random_types_values),
             numpy.average(nonrandom_types_values),
             scipy.stats.ttest_ind(random_types_values,
                                   nonrandom_types_values)[1]))

    text_output.append('Average number of mutated genes per random, non-randomly mutated strain: %f, %f, P-value (t-test, two-tailed): %g' % (
        numpy.average(random_gene_counts),
        numpy.average(nr_gene_counts),
        scipy.stats.ttest_ind(random_gene_counts, nr_gene_counts)[1]))

    random_go_processes = sql_interface.mutation_go_tags_by_doe(types_to_include=random_types)
    nr_go_processes = sql_interface.mutation_go_tags_by_doe(types_to_include=non_random_types)

    text_output.append('Average number of gene ontology processes changed per random, non-randomly mutated strain: %f, %f, P-value (t-test, two-tailed): %g' % (
        numpy.average(random_go_processes),
        numpy.average(nr_go_processes),
        scipy.stats.ttest_ind(random_go_processes, nr_go_processes)[1]))

    random_essential = sql_interface.gene_essentiality_by_doe(types_to_include=random_types)
    nr_essential = sql_interface.gene_essentiality_by_doe(types_to_include=non_random_types)

    text_output.append('Average number of essential genes modified per random, non-random strain: %f, %f, P-value (t-test, two-tailed): %g' \
          % (numpy.average(random_essential),
             numpy.average(nr_essential),
             scipy.stats.ttest_ind(random_essential, nr_essential)[1]))

    with open(os.path.join(constants.OUTPUT_DIR, 'Random vs. Non-random Analysis.txt'), 'w') as f:
        for line in text_output:
            f.write(line + '\n')


def phenotype_phenotype_interaction_analysis():

    _, ecoli_genome_size = sql_interface.get_distinct_entities('gene')
    association_btw_phenotypes(genome_size=ecoli_genome_size,
                               identifier='gene',
                               output_prefix='TagBased_',
                               ignore_threshold=4)

    _, ecoli_genome_size = sql_interface.get_distinct_entities('go')
    association_btw_phenotypes(genome_size=ecoli_genome_size,
                               identifier='go',
                               output_prefix='GO_',
                               ignore_threshold=4)

    _, ecoli_genome_size = sql_interface.get_distinct_entities('metabolite')
    association_btw_phenotypes(genome_size=ecoli_genome_size,
                               identifier='metabolite',
                               output_prefix='MSBased_',
                               ignore_threshold=4)


def error_analysis():

    mutations_to_process = {'aa_snps', 'nuc_snps'}
    gene_mutations = sql_interface.get_gene_mutation_tuples(filter_for_mutation_types=mutations_to_process)

    """

    gene_mutations = [{'strain': 'mg1655',
                       'mutation_type': 'nuc_snps',
                       'species_name': 'thrL',
                       'annotation': {'nuc_snps': [(189, 'A', 'A'),
                                                   (190, 'T', 'T'),
                                                   (191, 'G', 'G')]}}]
                                                   
    """

    genome_sequences = dict()
    for species in constants.SPECIES_LIST:
        genome_sequence = sql_interface.get_genome_sequence(query_species=species)
        genome_sequences[species] = genome_sequence

    nuc_error = 0
    nuc_total = 0
    nuc_backing_missing = 0

    aa_error = 0
    aa_total = 0
    aa_backing_missing = 0

    for result in gene_mutations:

        gene = result['species_name']

        strain = constants.get_strain_converter(result['strain'])

        if result['mutation_type'] not in mutations_to_process:
            continue

        # do sequence error analysis; the other types of data entry errors are hard to detect
        if 'nuc_snps' in result['mutation_type']:

            if strain not in genome_sequences:
                raise AssertionError('Genetic information missing for strain %s' % strain)

            # todo: add chromosome specifications to locations (indel, nucleotide, etc...)
            # correct specification
            if isinstance(genome_sequences[strain], list) and len(genome_sequences[strain]) > 1:
                raise AssertionError('Multi-chromosome case not handled for strain %s' % strain)

            gen_seq = genome_sequences[strain]

            for (location, original_base, altered_base) in result['annotation']['nuc_snps']:

                nuc_total += 1

                # incorrectly specified location
                if location is None or location > len(gen_seq):

                    nuc_backing_missing += 1
                    continue

                # base mismatch
                if gen_seq[location] != original_base:
                    nuc_error += 1

        if 'aa_snps' in result['mutation_type']:

            aa_seq = sql_interface.get_polypeptide_sequence(gene, query_species=strain)

            if aa_seq is None:
                aa_backing_missing += 1
                continue

            for (location, original_residue, altered_residue) in result['annotation']['aa_snps']:

                aa_total += 1

                if location is None or location > len(aa_seq):
                    aa_error += 1
                    continue

                if location != len(aa_seq) and aa_seq[location] != original_residue:
                    aa_error += 1

    with open(os.path.join(constants.OUTPUT_DIR, 'Error Results.txt'), 'w') as f:

        f.write('Nucleotide errors: %f' % (float(nuc_error + nuc_backing_missing) / float(nuc_total),) + '\n')
        f.write('AA errors: %f' % (float(aa_error + aa_backing_missing) / float(aa_total),) + '\n')


def generate_mutation_frequency_table():
    """
    
    Generates two tables:
    
    1. A list of every mutated gene in the Resistome, the associated phenotype, and mutations affecting it.
    2. A list of every mutated gene in the Resistome and the number of times it is mutated.
    
    :return: 
    """

    aggregate_table_data = sql_interface.get_aggregate_table_data()

    gene_counter = defaultdict(int)
    complete_mutation_report = []
    gene_counter_report = []

    associated_phenotypes = defaultdict(set)
    mg1655_name_dict = defaultdict(str)

    names = set()

    for result in aggregate_table_data:

        name = result['species_name'][0]
        mg_name = result['mg1655_name'][0]
        mut_types = result['mutation_type']
        annotations_details = result['annotation']
        phenotypes = result['phenotype']
        phenotype_types = result['phenotype_type']

        # count both to give users a feel for which species are actually mutated
        gene_counter[mg_name] += 1

        names.add(name)

        mg1655_name_dict[name] = mg_name

        seen_before = set()
        for phenotype, ptype in zip(phenotypes, phenotype_types):
            phenotype_type = 'Sensitive' if ptype == 'S' else 'Resistant'
            if (phenotype, phenotype_type) in seen_before:
                continue
            else:
                seen_before.add((phenotype, phenotype_type))

            associated_phenotypes[name].add((phenotype, phenotype_type))
            associated_phenotypes[mg_name].add((phenotype, phenotype_type))

        phenotypes = ';'.join([':'.join(x) for x in seen_before])
        annotation_text = []

        annotation_dict = {}
        for adict in annotations_details:
            for x in adict:
                annotation_dict[x] = adict[x]

        for mut_type in mut_types:
            annotation_text.append(mut_type + ':' + database_utils.standard_mutation_formatting(mut_type,
                                                                                                annotation_dict))
        complete_mutation_report.append((name, mg_name, phenotypes, ';'.join(mut_types), ';'.join(annotation_text)))

    for name in names:
        p_t_tuples = [x[0] + ':' + x[1] for x in associated_phenotypes[name]]

        gene_counter_report.append((name,
                                    gene_counter[mg1655_name_dict[name]],
                                    ';'.join(p_t_tuples)))

    with open(os.path.join(constants.OUTPUT_DIR, 'Complete Mutation Report.txt'), 'w') as f:

        f.write('\t'.join(['Name', 'MG1655 Accession', 'Phenotype', 'Mutation Type', 'Annotation']) + '\n')
        for x in complete_mutation_report:
            f.write('\t'.join(map(str, x)) + '\n')

    with open(os.path.join(constants.OUTPUT_DIR, 'Gene Count Report.txt'), 'w') as f:

        f.write('\t'.join(['Gene Name', 'Mutation Count', 'Associated Phenotypes']) + '\n')
        for x in gene_counter_report:
            f.write('\t'.join(map(str, x)) + '\n')


def cluster_statistics():

    gene_graph, _, _ = sql_interface.get_gene_gene_graph(edge_count_threshold=1)

    gene_to_cluster, cluster_membership = network_utils.compute_igraph_clustering(gene_graph, method='fastgreedy')

    common_gene_names = sql_interface.get_gene_common_names(gene_to_cluster.keys(), query_species='mg1655')

    with open(os.path.join(constants.OUTPUT_DIR, 'Cluster memberships.txt'), 'w') as f:

        for cluster in cluster_membership:

            if len(cluster_membership[cluster]) == 1:
                continue

            f.write('Cluster ID: %s\n' % str(cluster))
            output_array = []
            for gene in cluster_membership[cluster]:
                real_name = common_gene_names.get(gene, gene)
                output_array.append('%s (%s)' % (gene, real_name))
            f.write(', '.join(sorted(output_array)) + '\n')


def single_mutation_genotype_analysis():

    genotypes = sql_interface.get_mutant_genotypes()

    gene_phenotype_counter = defaultdict(int)

    multiple_counter = 0
    single_counter = 0

    for r in genotypes:

        if len(r['mg1655_names']) > 1:
            multiple_counter += 1
            continue

        single_counter += 1

        for gene in r['mg1655_names']:
            gene_phenotype_counter[gene] += len(r['phenotypes'])

    single_phenotype_association = sum([1 if gene_phenotype_counter[x] == 1 else 0 for x in gene_phenotype_counter])
    multiple_phenotype_association = sum([1 if gene_phenotype_counter[x] > 1 else 0 for x in gene_phenotype_counter])

    glib.histogram(list(gene_phenotype_counter.values()),
                   'Gene-Phenotype Associations (per gene)',
                   'Count',
                   os.path.join(constants.OUTPUT_DIR,
                                'Gene-Pheno Count Distribution.pdf'),
                   use_median_filter=True)


def run_analysis():

    # todo: distribution of mutation count
    # todo: phenotype distribution

    from resistome.machine_learning import recommendation_system

    error_analysis()

    vital_statistics()

    phenotype_phenotype_interaction_analysis()

    recommendation_system.generate_serialized_output()

    single_mutation_genotype_analysis()
    #
    # # set of regulators in MG1665
    # # list of go_tags related to regulation
    regulation_tags = ['GO:0016987', 'GO:0006351', 'GO:0006353']
    #
    # get list of genes that are regulators in the MG1655 dataset
    regulators = sql_interface.convert_features_to_genes(regulation_tags,
                                                        'go',
                                                         query_species='mg1655')

    map_genes_by_tag(regulators,
                     os.path.join(constants.OUTPUT_DIR,
                                  'Highly Connected Regulator HeatMap.pdf'),
                     tag_threshold=0,
                     aggregate_threshold=20,
                     binary=False)

    distribution_data = map_genes_by_tag(regulators,
                         os.path.join(constants.OUTPUT_DIR,
                                      'Completed Connected Regulator HeatMap.pdf'),
                         tag_threshold=0,
                         aggregate_threshold=0,
                         binary=False)

    glib.histogram(distribution_data.values(),
                   'Regulator-Phenotype Associations',
                   'Counts',
                   os.path.join(constants.OUTPUT_DIR,
                                'Distribution of Regulator-Phenotype counts.pdf'))

    generate_mutation_frequency_table()

    basic_information_figures()

    tag_distribution('Distribution of Resistome Tags.pdf', top=12)

    single_mutation_genotype_analysis()

    random_versus_designed_comparison()


if __name__ == '__main__':

    run_analysis()
