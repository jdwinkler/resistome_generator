from collections import defaultdict
import os
from resistome.sql import sql_interface
import numpy
from resistome.graphics import visualization as glib
from resistome.utils import go_enrichment
from resistome import constants


def get_unique_gene_residue_combinations():

    """
    
    Extracts tuples of (position, original aa, mutated aa) from the database while also returning the AA sequences
    corresponding to the mutated genes.
    
    :param cursor: 
    :return: 
    """

    mutation_data = sql_interface.get_gene_mutation_tuples(['aa_snps'])

    gene_mutation_tuples = set()

    for result in mutation_data:

        gene = result['species_name']
        mg_accession = result['mg1655_name']
        aa_snps = result['annotation']['aa_snps']
        strain = result['strain']

        aa_snp_set = set()

        for (position, oa, ma) in aa_snps:

            # these are stop codons and hard to handle in most analyses
            # usually assume to be inactivating if within coding sequence
            # removal of final stop has an unclear effect (also rare)
            if oa == '*' or ma == '*':
                continue

            aa_snp_set.add((position, oa, ma))

        gene_mutation_tuples.add((gene, mg_accession, frozenset(aa_snp_set), strain))

    sequence_cache = dict()
    filtered_tuples = []

    for (gene, mg_accession, mutations, strain) in gene_mutation_tuples:

        if gene in sequence_cache:
            aa_seq = sequence_cache[gene]
        else:
            aa_seq = sql_interface.get_polypeptide_sequence(gene, query_species=strain)

            # unknown genes
            if aa_seq is None:
                continue
            sequence_cache[gene] = aa_seq

        passes_qc = True

        for (position, oa, _) in mutations:

            # this means that the position is valid and the AA at position matches
            # that specified in the Resistome
            if position is None or position >= len(aa_seq) or aa_seq[position] != oa:
                passes_qc = False

        if passes_qc:
            filtered_tuples.append((gene, mg_accession, mutations, strain))

    return filtered_tuples, sequence_cache


def prepare_inps_files(filtered_tuples, seq_cache):

    """
    
    Prepares INPS files for predicting mutational effects.
    
    :param filtered_tuples: 
    :param seq_cache: 
    :return: 
    """

    sequence_mutations = defaultdict(set)

    for (gene, _, mutations, _) in filtered_tuples:

        sequence_mutations[gene].update(mutations)

    for gene in sequence_mutations:

        sequence = seq_cache[gene]
        mutations = []

        for (position, oa, ma) in sequence_mutations[gene]:
            mutations.append(oa + str(position+1) + ma)

        if len(mutations) == 0:
            continue

        with open(os.path.join(constants.ERROR_DIR, gene + '.fa'), 'w') as f:
            f.write('>' + gene + '\n')
            f.write(sequence + '\n')

        with open(os.path.join(constants.ERROR_DIR, gene + '.mut.list.txt'), 'w') as f:
            for mut in mutations:
                f.write(mut + '\n')


def analyze_inps_mutation_effects():

    """
    
    Summarizes the protein effects that are predicted to be associated with AA changes in the database
    by INPS.
    
    :param cursor: 
    :return: 
    """

    protein_stability_data = sql_interface.extract_protein_stability_data('INPS')

    stability_values = []

    stabilized_proteins = []

    destabilized_count = 0
    stabilized_count = 0

    for result in protein_stability_data:

        stability_values.append(result['score'])

        if result['score'] < 0.0:
            stabilized_count += 1
        else:
            destabilized_count += 1

        # stabilizing ddG change
        if result['score'] <= -1.0:
            stabilized_proteins.append(result)

    return stability_values, stabilized_count, destabilized_count


def filter_aa_tuples_for_snap2(filtered_tuples, number_of_top_genes):

    """
    
    Identifies the genes with the most AA changes in the resistome dataset. Serves to decouple the SNAP2 heatmap
    method from extracting genes of interest.
    
    :param filtered_tuples: 
    :param number_of_top_genes: 
    :return: 
    """

    count_dict = defaultdict(int)

    for (gene, mg_accession, aa_mutations, host_strain) in filtered_tuples:

        count_dict[mg_accession] += 1

    count_list = sorted([(x, y) for x, y in count_dict.items()], key=lambda x: x[1], reverse=True)

    genes_to_process = set([x[0] for x in count_list[0:number_of_top_genes]])

    return genes_to_process


def generate_snap2_frequently_hit_heatmaps(filtered_tuples,
                                           genes_to_process,
                                           sequence_dict,
                                           force_identical_sequences=False,
                                           force_specific_strain=False,
                                           target_strain='mg1655',
                                           skip_errors=False):

    """
    
    Identifies the most heavily mutated (by residue changes) proteins in the database and uses SNAP2 to:
    
    figure out how much of the mutational space has been sampled
    the predicted effects of these mutations
    a matrix for visualizing these in heatmap form
    
    Outputs heat maps with marked sites indicating mutations.
    
    Note that force identical sequences may result in no output for a given call; any non-identical genes are ignore
    since you cannot mix their AA mutation data easily. You can specify a specific strain (mg1655, rel606, or w)
    if you want to limit the analysis and avoid this problem, but all data from other strains will be discarded.
    
    You can skip errors to effectively ignore sequence differences (any residues that differ and have mutations will
    skipped).
    
    :param filtered_tuples:
    :param genes_to_process:
    :param sequence_dict:
    :param force_identical_sequences:
    :param force_specific_strain:
    :param target_strain:
    :return: 
    """

    if target_strain not in constants.SPECIES_LIST:
        raise AssertionError('Unknown strain')

    aa_list = set(constants.AMINO_ACIDS) - {'*'}

    aa_list = sorted(aa_list, reverse=True)

    filtered_tuples = list(filter(lambda x: x[1] in set(genes_to_process), filtered_tuples))

    for gene in genes_to_process:

        heatmap_matrix = []
        collated_tuples = []
        processed_set = set()
        rectangles_to_highlight = []

        temp_filtered_tuples = list(filter(lambda x: x[1] == gene, filtered_tuples))
        species_key = set([x[0] for x in temp_filtered_tuples])

        identical = True
        for sk in species_key:
            for sy in species_key:
                if sk == sy:
                    continue
                if sequence_dict[sk] != sequence_dict[sy]:
                    identical = False

        if not identical and force_identical_sequences:
            continue
        else:
            sequence_identifier = species_key.pop()

        for (species_accession, mg_accession, aa_mutations, strain) in temp_filtered_tuples:

            if force_specific_strain and target_strain != strain:
                continue

            if force_specific_strain:
                gene = species_accession

            collated_tuples.extend(aa_mutations)

        # possible to have no data to process for this gene
        if len(collated_tuples) == 0:
            continue

        tuples_to_query = []
        for aa_row, row_position in zip(aa_list, range(0, len(aa_list))):

            for aa_col, protein_position in zip(sequence_dict[sequence_identifier],
                                        range(0, len(sequence_dict[sequence_identifier]))):

                if (protein_position, aa_col, aa_row) in collated_tuples:

                    processed_set.add((protein_position, aa_col, aa_row))
                    tuples_to_query.append((protein_position, aa_col, aa_row))

                elif aa_row != aa_col:
                    tuples_to_query.append((protein_position, aa_col, aa_row))

        query_result = sql_interface.get_aa_stability_position_value(gene=gene,
                                                                     target_strain='mg1655',
                                                                     positions_wt_mut_tuples=tuples_to_query,
                                                                     method='SNAP2',
                                                                     skip_errors=skip_errors)

        for aa_row, row_position in zip(aa_list, range(0, len(aa_list))):

            aa_row_data = []
            for aa_col, protein_position in zip(sequence_dict[sequence_identifier],
                                                range(0, len(sequence_dict[sequence_identifier]))):

                if (protein_position, aa_col, aa_row) in query_result:
                    aa_row_data.append(query_result[(protein_position, aa_col, aa_row)])
                    if (protein_position, aa_col, aa_row) in collated_tuples:
                        rectangles_to_highlight.append((row_position, protein_position,
                                                        query_result[(protein_position, aa_col, aa_row)], 'yellow'))
                else:
                    aa_row_data.append(0)

            heatmap_matrix.append(aa_row_data)

        # all tuples should be correct due to validation in get_unique_gene_residue
        # EXCEPT if your gene is missing from the SNAP2 database

        # assert len(set(collated_tuples) - processed_set) == 0, 'Did not process all tuples!'
        
        # mutations are highlighted with * (asterixs)
        # includes color bar
        # positive values predict effects, negative values neutrality
        glib.generate_heatmap_asymmetric(matrix=heatmap_matrix,
                                         x_names=sequence_dict[sequence_identifier],
                                         y_names=aa_list,
                                         ylabel='SNAP2 Fitness Estimate',
                                         filename=os.path.join(constants.OUTPUT_DIR,
                                                               gene + '_SNAP2_replacement_map.pdf'),
                                         plot_color_bar=True,
                                         dim_override=((len(sequence_dict[sequence_identifier])/6,
                                                        len(aa_list)/5)),
                                         cmap='seismic',
                                         patches=rectangles_to_highlight,
                                         vmax=1.0,
                                         vmin=-1.0)


def protein_domain_analysis(domains_to_ignore=None):

    """
    
    Calculates the statistical enrichment of protein domains that are affected by AA changes in the database. Returns
    the number of times these regions are mutated (rv), the p-value of their enrichment (region_pvalues), genes
    missing data, and the number of mutations found.
    
    :param domains_to_ignore:
    :return: 
    """

    domains_to_ignore = set()

    mutation_data = sql_interface.get_gene_mutation_tuples(['aa_snps'])
    missing_data = 0
    region_counter = defaultdict(int)

    for result in mutation_data:

        gene_name = result['mg1655_name']
        strain = result['strain']

        annotation = result['annotation']['aa_snps']

        for (location, _, _) in annotation:

            # returns domains that contain this residue
            domains = sql_interface.get_uniprot_location_data(strain, gene_name, location)
            domains = list(filter(lambda x: x['region'] not in domains_to_ignore, domains))

            if len(domains) == 0:
                missing_data += 1
            else:
                for x in domains:
                    region_counter[x['region']] += 1

    region_values = []

    for region in region_counter:

        region_values.append((region, region_counter[region]))

    rv = sorted(region_values, key=lambda x: x[1], reverse=True)

    # fisher's exact test to evaluate hits to a given protein domain
    e_calculator = go_enrichment.EnrichmentCalculator(sql_interface.get_gene_label_relationships('uniprot'))

    region_pvalues = e_calculator.compute_enrichment_specify_counts(region_counter,
                                                                    constants.P_VALUE_THRESHOLD)

    return rv, region_pvalues, missing_data, len(mutation_data)


def run_analysis():

    stability_values, stabilized_count, destabilized_count = analyze_inps_mutation_effects()

    output_text = []

    output_text.append('Median ddG: %g' % (numpy.median(stability_values),))
    output_text.append('Mean ddG: %g' % (numpy.mean(stability_values),))
    output_text.append('Ratio of stabilized to destabilized proteins: %g'
                       % (float(stabilized_count) / float(destabilized_count),))

    filtered_tuples, seq_dict = get_unique_gene_residue_combinations()
    genes_of_interest = filter_aa_tuples_for_snap2(filtered_tuples, number_of_top_genes=15)

    # can also specify specific genes of interest, top 10 here just as an example essentially
    # also quite slow
    generate_snap2_frequently_hit_heatmaps(filtered_tuples,
                                           genes_of_interest,
                                           seq_dict,
                                           skip_errors=True)

    glib.histogram(stability_values, 'ddG','Count', os.path.join(constants.OUTPUT_DIR, 'Stability Histogram.pdf'))

    domains_to_ignore = {'Chain'}

    region_values, region_p, missing, results = protein_domain_analysis(domains_to_ignore=domains_to_ignore)

    normalizing_sum = float(sum([x[1] for x in region_values[1:]]))

    output_text.append('Fraction missing backing protein structural data: %g' % (float(missing)/float(results)))

    glib.bargraph([x[0] for x in region_values[1:]],
                  [float(x[1])/normalizing_sum for x in region_values[1:]],
                         'Protein Domain Type',
                         'Count',
                  os.path.join(constants.OUTPUT_DIR,
                                      'Protein Domain Mutagenesis.pdf'),
                  rotation='vertical')

    output_text.append('\nRegion\tEnrichment P-value\tOdds Ratio')

    for (region, p, odds) in region_p:
        output_text.append('\t'.join((region, str(p), str(odds))))

    with open(os.path.join(constants.OUTPUT_DIR,
                           'Protein Domain Statistical Analysis.txt'), 'w') as f:
        for line in output_text:
            f.write(line + '\n')


if __name__ == '__main__':

    run_analysis()