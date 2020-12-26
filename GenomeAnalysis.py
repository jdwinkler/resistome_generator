from resistome import constants
import numpy
from resistome.graphics import visualization as glib
from collections import defaultdict
import scipy.stats
import os
from resistome.sql import sql_interface


def hamming_distance(str1, str2):

    """
    
    Calculates the number of different characters between two strings of the same length.
    
    Throws an AssertionError if you pass unequal length strings.
    
    :param str1: 
    :param str2: 
    :return: 
    """

    distance = 0
    differences = []
    if len(str1) != len(str2):
        raise AssertionError('unequal length in strings')

    # calculate hamming distance
    for base_x, base_y in zip(str1, str2):
        if base_x != base_y:
            distance += 1
            differences.append((base_x, base_y))
        else:
            # identical bases
            differences.append(('*', '*'))

    return distance, differences


def gatc_analysis(window=200):

    """
    
    Examines the distributions of GATC sites around the genome and the associated SNP/AA replacement rates
    in the general vicinity. Note that this analysis is very sensitive to the window size considered, and hence
    not very robust/reliable. 
    
    Window refers to the amount of flanking sequence to check for GATC sites around a given mutation.
    
    :param window: 
    :return: 
    """

    location_data = sql_interface.get_gene_location_data(query_species='mg1655')
    genome = sql_interface.get_genome_sequence(query_species='mg1655')

    gene_counter = defaultdict(int)

    gene_mutation_tuples = sql_interface.get_gene_mutation_tuples(filter_for_mutation_types=['indel',
                                                                                            'nuc_snps',
                                                                                            'aa_snps'])
    for result in gene_mutation_tuples:

        gene_counter[result['mg1655_name']] += 1

    mutation_gatc_correlation = []

    for gene in gene_counter:

        if gene not in location_data or gene_counter[gene] <= 1:
            continue

        start = location_data[gene]['start']
        end = location_data[gene]['stop']

        # count GATC sites in flanking sequences before coding
        before_gene_gatc = genome[start - window:start].count('GATC') + genome[start - window:start].count('CTAG') \
            if start - window > 0 else 0

        # same, but within coding region
        coding_region_gatc = genome[start:end].count('GATC') + genome[start:end].count('CTAG')

        # after
        after_gene_gatc = genome[end:end + window].count('GATC') + genome[end:end + window].count('CTAG') \
            if end + window < len(genome) else 0

        gatc_count = float(before_gene_gatc + coding_region_gatc + after_gene_gatc)

        mutation_gatc_correlation.append((gene_counter[gene], gatc_count))

    return mutation_gatc_correlation


def snp_spectrum_analysis():

    """
    
    Generates a matrix of data plus transition counts (N => N' : int) describing nucleotide substitution patterns
    in the Resistome.
    
    Note that the location data for the SNPs does not seem accurate, but it is assumed that the actual base changes
    are. Also includes AA => AA' transitions since nearly all of those are characterized by SNPs too.
    
    :return: 
    """

    import itertools

    code, aa_to_codon = sql_interface.get_genetic_code()
    # extract and munge all SNP data
    annotation_data = sql_interface.get_gene_mutation_tuples(['nuc_snps', 'aa_snps'])

    aa_mutations = []
    nuc_mutations = []

    for result in annotation_data:
        aa_mutations.extend(result['annotation'].get('aa_snps', []))
        nuc_mutations.extend(result['annotation'].get('nuc_snps', []))

    # accounts for both aa and nuc changes
    snp_transition_counts = defaultdict(dict)
    amino_specific_counts = defaultdict(dict)
    nucleotide_specific_counts = defaultdict(dict)

    for (pos, aa_x, aa_y) in aa_mutations:

        # more standard nomenclature
        if aa_x == '*':
            aa_x = 'X'
        if aa_y == '*':
            aa_y = 'X'
            
        x_codons = aa_to_codon[aa_x]
        y_codons = aa_to_codon[aa_y]
        cp = itertools.product(x_codons, y_codons)
        min_dist_results = []

        for (x, y) in cp:
            d, differences = hamming_distance(x, y)
            # AA => AA' changes are only included if they can be caused by SNPs alone (most ~98% are)
            if d == 1:
                min_dist_results.append(differences)

        for result in min_dist_results:
            for (source, sink) in result:
                if source == '*':
                    continue
                if sink in snp_transition_counts[source]:
                    snp_transition_counts[source][sink] += 1
                    amino_specific_counts[source][sink] += 1
                else:
                    snp_transition_counts[source][sink] = 1
                    amino_specific_counts[source][sink] = 1

    for (pos, original_base, mutated_base) in nuc_mutations:

        if mutated_base in snp_transition_counts[original_base]:
            snp_transition_counts[original_base][mutated_base] += 1
        else:
            snp_transition_counts[original_base][mutated_base] = 1

        if mutated_base in nucleotide_specific_counts[original_base]:
            nucleotide_specific_counts[original_base][mutated_base] += 1
        else:
            nucleotide_specific_counts[original_base][mutated_base] = 1

    keys = ['A', 'T', 'G', 'C']

    matrix = []

    normalizing_sum = 0
    for x in keys:
        for y in keys:
            normalizing_sum += snp_transition_counts[x][y] if y in snp_transition_counts[x] else 0

    normalizing_sum = float(normalizing_sum)

    for base_x in keys:
        temp = []
        for base_y in keys:
            value = snp_transition_counts[base_x][base_y] if base_y in snp_transition_counts[base_x] else 0
            temp.append(float(value))
        matrix.append([x/normalizing_sum for x in temp])

    glib.generate_heatmap(matrix,
                          keys,
                          'N->N\' Conversion',
                          os.path.join(constants.OUTPUT_DIR, 'Mutation N to N Heatmap.pdf'),
                          cmap='viridis')

    bar_graph_data = []

    for x in keys:
        for y in keys:
            if x == y:
                continue
            bar_graph_data.append((x + '>' +y, snp_transition_counts[x].get(y, 0)))

    bar_graph_data = sorted(bar_graph_data, key=lambda x: x[0])
    mutation_total = sum([y[1] for y in bar_graph_data])

    glib.bargraph([x[0] for x in bar_graph_data],
                  [float(y[1])/float(mutation_total) for y in bar_graph_data],
                  'Nucleotide Transition (N>N\')',
                  'Fraction of Observed SNPs',
                  os.path.join(constants.OUTPUT_DIR,
                               'Mutation N to N Bargraph.pdf'),
                  rotation='vertical')

    bar_graph_data = []

    for x in keys:
        for y in keys:
            if x == y:
                continue
            bar_graph_data.append((x + '>' +y, amino_specific_counts[x].get(y, 0)))

    bar_graph_data = sorted(bar_graph_data, key=lambda x: x[0])
    aa_mutation_total = sum([y[1] for y in bar_graph_data])

    glib.bargraph([x[0] for x in bar_graph_data],
                  [float(y[1])/float(aa_mutation_total) for y in bar_graph_data],
                  'Nucleotide Transition (N>N\')',
                  'Fraction of Observed SNPs',
                  os.path.join(constants.OUTPUT_DIR,
                               'Mutation N to N Bargraph (AA only).pdf'),
                  rotation='vertical')

    bar_graph_data = []

    for x in keys:
        for y in keys:
            if x == y:
                continue
            bar_graph_data.append((x + '>' +y, nucleotide_specific_counts[x].get(y, 0)))

    bar_graph_data = sorted(bar_graph_data, key=lambda x: x[0])
    snp_mutation_total = sum([y[1] for y in bar_graph_data])

    glib.bargraph([x[0] for x in bar_graph_data],
                  [float(y[1])/float(snp_mutation_total) for y in bar_graph_data],
                  'Nucleotide Transition (N>N\')',
                  'Fraction of Observed SNPs',
                  os.path.join(constants.OUTPUT_DIR,
                               'Mutation N to N Bargraph (Explicit SNPs only).pdf'),
                  rotation='vertical')

    bar_graph_data = []

    for x in keys:
        for y in keys:
            if x == y:
                continue

            aa_value = float(amino_specific_counts[x].get(y, 0))/float(aa_mutation_total)
            nuc_value = float(nucleotide_specific_counts[x].get(y, 1))/float(snp_mutation_total)

            bar_graph_data.append((x + '>' +y, float(aa_value)/float(nuc_value)))

    bar_graph_data = sorted(bar_graph_data, key=lambda x: x[0])

    glib.bargraph([x[0] for x in bar_graph_data],
                  [y[1] for y in bar_graph_data],
                  'Nucleotide Transition (N>N\')',
                  'Ratio of AA Inferred to Explicit SNPs',
                  os.path.join(constants.OUTPUT_DIR,
                               'Mutation N to N Bargraph (Ratio).pdf'),
                  rotation='vertical')

    return matrix, snp_transition_counts, amino_specific_counts, nucleotide_specific_counts


def context_analyzer(context_sequence, context_length):

    """
    
    Analyzes the homopolymeric content of context sequence (mainly for indel analysis).
    
    :param context_sequence: 
    :param context_length: 
    :return: 
    """

    def contains_min_homopolymer(seq, run=3):

        """
        
        Detects runs of length 'run' in seq.
        
        :param seq: 
        :param run: 
        :return: 
        """

        return 'A' * run in seq \
               or 'T' * run in seq \
               or 'C' * run in seq \
               or 'G' * run in seq

    distance_from_indel = defaultdict(list)

    triplet_counter = defaultdict(int)

    # ^ is the indel start site
    min_len = min([len(x[0].replace('^', '')) for x in context_sequence])

    run_containing_count = 0
    triplets = 0
    coding_sequence_count = 0
    genomic_indel_count = 0

    for (sequence, codons) in context_sequence:

        caret_index = sequence.find('^')
        indel_location = caret_index + 1
        surrounding_sequence = sequence[indel_location - 3:indel_location + 4].replace('^', '')

        if contains_min_homopolymer(surrounding_sequence):
            run_containing_count += 1

        # looks at codons that are affected by the indel (potentially)

        if len(codons) == 0:

            fake_codon_start = sequence[indel_location-4:indel_location-1]
            fake_codon_middle = sequence[indel_location:indel_location+3]

            end_pos = 3

            fake_codon_end = sequence[indel_location+end_pos:indel_location+end_pos+3]

            if '^' in fake_codon_middle:
                end_pos = 4
                fake_codon_middle = sequence[indel_location:indel_location+4].replace('^', '')
                fake_codon_end = sequence[indel_location+end_pos:indel_location+end_pos+3]

            if '^' in fake_codon_end:

                if fake_codon_end[0] == '^':
                    fake_codon_end = sequence[indel_location+end_pos+1:indel_location+end_pos+4]
                if fake_codon_end[1] == '^' or fake_codon_end[2] == '^':
                    fake_codon_end = sequence[indel_location+end_pos:indel_location+end_pos+4].replace('^', '')

            genomic_indel_count += 1
            coding_sequence_count -= 1

            codons = [fake_codon_start, fake_codon_middle, fake_codon_end]

        if len(codons) >= 2:
            coding_sequence_count += 1
            triplets += 1
            triplet_counter[codons[0] + '\n' + codons[1]] += 1
            if len(codons) == 3:
                triplet_counter[codons[1] + '\n' + codons[2]] += 1

        # remove caret
        sequence = sequence.replace('^', '')

        for i in range(1, min_len - 1):

            if sequence[i] == sequence[i - 1]:
                distance_from_indel[i].append(1)
            else:
                distance_from_indel[i].append(0)

    indel_dist_fraction = []

    # distance away from homopolymer/insert site
    for key in distance_from_indel:
        indel_dist_fraction.append((key - 1 - context_length, numpy.mean(distance_from_indel[key])))

    indel_dist_fraction = sorted(indel_dist_fraction, key=lambda x: x[0])

    triplet_fractions = []

    sum_triplets = float(sum(triplet_counter.values()))

    # fraction of triplets that are found near/with indels
    for key in triplet_counter:
        triplet_fractions.append((key, float(triplet_counter[key]) / sum_triplets))

    descriptive_stats_dict = {'genomic_indels': genomic_indel_count,
                              'total_indels': len(context_sequence),
                              'run_containing_indels': float(run_containing_count) / float(len(context_sequence))}

    return triplet_fractions, indel_dist_fraction, descriptive_stats_dict


def get_flanking_range(pos, size, genome_len, flank_len):

    """
    
    Extracts the flanking sequences around position with the specified size.
    
    Note that this method will throw an assertion error if a flank around the genome zero point is requested (e.g.
    if thrL is mutated, for example). Can be updated to fix this if needed.
    
    :param pos: 
    :param size: 
    :param genome_len: 
    :param flank_len: 
    :return: 
    """

    if pos - flank_len > 0:
        fstart = pos - flank_len
    else:
        raise AssertionError('Need to implement wrap-around mechanics')

    if pos + size + flank_len < genome_len:
        fend = pos + size + flank_len
    else:
        raise AssertionError('Need to implement wrap-around mechanics')

    return fstart, fend


def get_codon_context(genome_sequence,
                      gene_start,
                      gene_stop,
                      position_affected):

    """
    
    Extracts codons containing and surrounding a given mutation (position_affected).
    
    If the codon is at position 0 (start), only the first and second codons are returned
    If the codon is at position N (end), only the second to last and last codons are returned
    otherwise the previous codon, affected codon, and the next codon are returned
    
    :param genome_sequence: 
    :param gene_start: 
    :param gene_stop: 
    :param position_affected: 
    :return: 
    """

    gene_length = gene_stop - gene_start + 1

    if gene_length % 3 != 0:
        # usually RNA of some sort
        return []

    if position_affected < gene_start:
        # intergenic mutation, no codons to return
        return []

    if position_affected > gene_stop:
        # same as above just after the gene
        return []

    # split into codon
    codons = [genome_sequence[i:i + 3] for i in range(gene_start, gene_stop + 1, 3)]
    codon_translation_dict = {}

    counter = 0

    # translate position => codon
    for codon in codons:

        codon_translation_dict[counter * 3 + 0 + gene_start] = counter
        codon_translation_dict[counter * 3 + 1 + gene_start] = counter
        codon_translation_dict[counter * 3 + 2 + gene_start] = counter

        counter += 1

    # stop codon
    if position_affected not in codon_translation_dict:
        return []

    codon_affected = codon_translation_dict[position_affected]

    # handle edge (literally) cases
    if codon_affected == 0:
        return [codons[0], codons[1]]
    if codon_affected == len(codons) - 1:
        return [codons[-2], codons[-1]]

    return [codons[codon_affected-1], codons[codon_affected], codons[codon_affected+1]]


def get_genome_context(genome_sequence,
                       gene_descriptor_dict,
                       gene,
                       location,
                       context_length=20,
                       max_indel_size=5):

    """
    
    Extract the genome (sequence) context surrounding a given indel. Mainly for GATC analysis.
    
    :param genome_sequence: 
    :param gene_descriptor_dict: 
    :param gene: 
    :param location: 
    :param context_length: 
    :param max_indel_size: 
    :return: 
    """

    allowed_types = {'absolute', 'absolute_inferred'}

    sequences = []

    (position, size, type_of_indel) = location[0], location[1], location[2]

    if position is None:
        return []

    if type_of_indel not in allowed_types:
        return []

    # descriptive location
    if isinstance(position, str):
        return []

    # fails to pass size restriction
    if max_indel_size is not None and abs(size) > max_indel_size:
        return []

    # find position in genome, take context length flanking sequences.
    if position > len(genome_sequence) or position < 0:
        return []

    if gene not in gene_descriptor_dict:
        return []

    size = abs(size)

    start = gene_descriptor_dict[gene]['start']
    stop = gene_descriptor_dict[gene]['stop']

    codons = get_codon_context(genome_sequence,
                               start,
                               stop,
                               position)

    (fstart, fend) = get_flanking_range(position,
                                        size,
                                        len(genome_sequence),
                                        context_length)

    # '^' is the site of the indel
    seq = list(genome_sequence[fstart:position] + '^'
               + genome_sequence[position:position + size] + '^'
               + genome_sequence[position + size:fend])
    sequences.append((str(''.join(seq)), codons))

    return sequences


def get_gene_indel_locations():

    """
    
    Does some basic filtering to extract usable indel data for analysis; some indel data is a bit wonky
    in that some datasets do not include proper locations for the indel.
    
    :return: 
    """

    data_store = []

    annotation_data = sql_interface.get_gene_mutation_tuples(['indel'])

    seen = set()

    for result in annotation_data:

        data = result['annotation']['indel']

        gene_name = result['mg1655_name'].upper()

        # makes sure the location is an integer + an absolute location (not descriptive)
        data = [(gene_name, x) for x in data
                if (x[2] == 'absolute' or x[2] == 'absolute_inferred') and isinstance(x[1], int)]

        data_store.extend(data)

    return data_store


def extract_indel_location_information(gene_indel_data_store, context_length=9):

    """
    
    Exracts the genome context for indel data (limited to small indels to avoid grabbing large deletions).
    
    :param gene_indel_data_store: 
    :param context_length: 
    :return: 
    """

    genome_sequence = sql_interface.get_genome_sequence(query_species='mg1655')
    genome_data = sql_interface.get_gene_location_data(query_species='mg1655')

    contexts = []

    zero_count = 0

    for (gene_name, locations) in gene_indel_data_store:

        try:

            sequence_context = get_genome_context(genome_sequence,
                                                  genome_data,
                                                  gene_name,
                                                  locations,
                                                  context_length=context_length + 1,
                                                  max_indel_size=None)

            if len(sequence_context) == 0:
                zero_count += 1

            contexts.extend(sequence_context)
        except AssertionError:
            pass


    return contexts


def compute_average_mutational_distance():

    """
    
    Calculates the average, median, and standard deviation of distance between mutations in each Resistome genotype.
    Also outputs a histogram (../output/Mutational Separation within Genotypes.pdf) of the data.
    
    :return: 
    """

    import itertools

    # select all genotypes with more than 1 mutation
    # this is actually a pretty crude way of doing this
    # only focuses on genes rather than exact mutation locations
    # get enough though?

    genotypes = sql_interface.get_mutant_genotypes()
    gene_locations = sql_interface.get_gene_location_data(query_species='mg1655')

    genome_length = float(len(sql_interface.get_genome_sequence(query_species='mg1655')))

    genotypes = list(filter(lambda x: len(x['species_names']) > 1, genotypes))

    mutation_pairwise_distances = []

    for genotype in genotypes:

        genes = genotype['mg1655_names']
        genes = filter(lambda x: x in gene_locations, genes)
        start_list = [gene_locations[x]['start'] for x in genes]

        # cartesian product of every start location
        start_product = itertools.product(start_list, start_list)
        observed_previously = set()

        for (s1, s2) in start_product:

            if (s1, s2) in observed_previously or s1 == s2:
                continue

            # TODO needs memory optimization.
            mutation_pairwise_distances.append(abs(s1 - s2)/1000.0)
            observed_previously.add((s1, s2))
            observed_previously.add((s2, s1))

    avg_distance = numpy.mean(mutation_pairwise_distances)
    median_distance = numpy.median(mutation_pairwise_distances)
    std_distance = numpy.std(mutation_pairwise_distances)

    glib.histogram(mutation_pairwise_distances,
                   'Mutation Separation (kilobases)',
                   'Counts',
                   os.path.join(constants.OUTPUT_DIR,
                                'Mutational Separation within Genotypes.pdf'))

    return avg_distance, median_distance, std_distance


def gene_hit_analysis(mutated_gene_names, genes_in_target_genome):

    """
    
    Screens the Resistome for genes not associated with any phenotype.
    
    :param mutated_gene_names: set of mutated gene names
    :param genes_in_target_genome: set of gene names
    :return: 
    """

    from resistome.utils.go_enrichment import EnrichmentCalculator

    calculator = EnrichmentCalculator(sql_interface.get_gene_label_relationships('go'))
    unmutated_genes = set(genes_in_target_genome) - set(mutated_gene_names)
    enriched_go_data = calculator.compute_enrichment(unmutated_genes, p_value_filter=constants.P_VALUE_THRESHOLD)

    return unmutated_genes, enriched_go_data


def run_analysis():

    """
    
    Driver method.
    
    :return: 
    """

    # avg_mut_dist, med_mut_dist, std_mut_dist = compute_average_mutational_distance()

    output_file = []

    gene_data_array = get_gene_indel_locations()

    contexts = extract_indel_location_information(gene_data_array, context_length=9)

    triplet_fractions, indel_dist_fraction, stats = context_analyzer(contexts, context_length=9)
    triplet_fractions = sorted(triplet_fractions, key=lambda x: x[1], reverse=True)[0:15]

    glib.bargraph([x[0] for x in indel_dist_fraction],
                  [y[1] for y in indel_dist_fraction],
                  'Distance from indel (bp)',
                  'P(Si == Si-1)',
                  os.path.join(constants.OUTPUT_DIR,
                               'Indel Transition Probabilities.pdf'))

    glib.bargraph([x[0] for x in triplet_fractions],
                  [y[1] for y in triplet_fractions],
                  'Indel Sequence Context',
                  'Frequency Near Indel [Adjacent or Affected Codon]',
                  os.path.join(constants.OUTPUT_DIR, 'Sequence Indel Contexts.pdf'))

    output_file.append('Fraction of indels containing at least trinucleotide run: %g' % stats['run_containing_indels'])
    output_file.append('Number of genomic indels: %i' % stats['genomic_indels'])
    output_file.append('Number of total indels: %i' % stats['total_indels'])

    # output_file.append('Mutational separation stats: %f  avg, %f median, %f, std dev' % (avg_mut_dist,
    #                                                                                      med_mut_dist,
    #                                                                                      std_mut_dist))

    window = 2000
    mutation_gatc_correlation = gatc_analysis(window=window)

    (d, p) = scipy.stats.spearmanr([x[0] for x in mutation_gatc_correlation],
                                   [x[1] for x in mutation_gatc_correlation])

    output_file.append('P-value for mutation frequency-gatc correlation: %g' % p)
    glib.scatter([numpy.log2(x[0]) for x in mutation_gatc_correlation],
                 [x[1] for x in mutation_gatc_correlation],
                 'Gene Hits',
                 'GATC Abundance in %i bp window' % window,
                 os.path.join(constants.OUTPUT_DIR,
                              'GATC Abundance vs Gene Mutation Frequency.pdf'))

    snp_spectrum_analysis()

    with open(os.path.join(constants.OUTPUT_DIR, 'Genome Analysis Report.txt'), 'w') as f:
        for line in output_file:
            f.write(line+'\n')

    gene_mutation_tuples = sql_interface.get_gene_mutation_tuples()
    mutated_gene_set = set([x['mg1655_name'] for x in gene_mutation_tuples])
    complete_gene_set = sql_interface.get_gene_identifiers()

    unmutated_genes, go_pval_tuples = gene_hit_analysis(mutated_gene_set, complete_gene_set)

    go_tag_name = sql_interface.get_go_terms([x[0] for x in go_pval_tuples])
    gene_names = sql_interface.get_gene_common_names(unmutated_genes)

    with open(os.path.join(constants.OUTPUT_DIR, 'Analysis of Non-Mutated Genes.txt'), 'w') as f:

        f.write('Genes not being mutated in the MG1655 genome:\n')
        f.write('Accession\tGene Name\n')

        for gene in unmutated_genes:
            f.write('\t'.join([gene, gene_names[gene]]) + '\n')

        f.write('\nGene ontology enrichment\n')
        f.write('Tag\tName\tP-value\n')

        for tag, p, fold_enrichment in go_pval_tuples:
            f.write('\t'.join([tag, go_tag_name[tag], str(p)]) + '\n')

if __name__ == '__main__':

    run_analysis()
