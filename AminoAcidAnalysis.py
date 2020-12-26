from resistome import constants
from resistome.graphics import visualization as glib
import os
import itertools
from collections import defaultdict
import numpy
import scipy.stats
from resistome.sql import sql_interface

__author__ = 'jdwinkler'


def compute_gene_fractions(gene_count):

    """
    
    Turns a dict of (AA_1 : dict (AA_2: list)) into a frequency dict by dividing the counts by the total
    number of genes affected in the dict.
    
    :param gene_count: 
    :return: 
    """

    gene_fraction = defaultdict(dict)

    all_genes = set()
    for aa_x in gene_count:
        for aa_y in gene_count:

            if aa_y in gene_count[aa_x]:
                all_genes.update(gene_count[aa_x][aa_y])

    norm = float(len(all_genes))

    for aa_x in gene_count:
        for aa_y in gene_count:
            if aa_y in gene_count[aa_x]:
                gene_fraction[aa_x][aa_y] = float(len(gene_count[aa_x][aa_y])) / norm
            else:
                gene_fraction[aa_x][aa_y] = 0

    return gene_fraction


def build_aa_transition_dict(aa_tuples, aa_distance_dict):

    """
    
    Computes the following information from:
    
    aa_tuples: list of (original amino acid, mutated amino acid, position)
    aa_distance_dict: minimum number of SNPs required to convert AA_1 to AA_2
    
    to:
    
    The number of transitions from a given amino acid
    The number of transitions from a given amino acid to another amino acid
    The distance of each transition from a given AA
    
    :param aa_tuples: 
    :param aa_distance_dict: 
    :return: 
    """

    paired_count_dict = defaultdict(dict)
    aa_transition_weighted = defaultdict(dict)

    overall_trans_count = defaultdict(int)

    for (oa, ma, pos) in aa_tuples:

        overall_trans_count[oa] += 1

        if ma in paired_count_dict[oa]:
            paired_count_dict[oa][ma] += 1
        else:
            paired_count_dict[oa][ma] = 1

        if aa_distance_dict[oa][ma] in aa_transition_weighted[oa]:
            aa_transition_weighted[oa][aa_distance_dict[oa][ma]] += 1
        else:
            aa_transition_weighted[oa][aa_distance_dict[oa][ma]] = 1

    return overall_trans_count, paired_count_dict, aa_transition_weighted


def generate_aa_stacked_figure(aa_order, total_aa_counter, aa_transition_weighted):

    """
    
    Generates a stacked bar figure indicating the fraction of AA transitions that are 1, 2, and 3 nucleotides
    away from a given AA. Entries that indicate a AA => AA transition (same AA) are ignored, though reported
    occasionally in the resistome.
    
    Outputs file to ../outputs/Estimated Nucleotide Transitions.pdf
    
    :param aa_order: 
    :param total_aa_counter: 
    :param aa_transition_weighted: 
    :return: 
    """

    datadict_array = []
    replacement_order = [1, 2, 3]

    for aa in aa_order:

        transitions = aa_transition_weighted[aa]

        for key in replacement_order:
            if key not in transitions:
                transitions[key] = 0

        value_sum = float(sum(transitions.values()))

        output_dict = {}
        for k in transitions:
            if k > 0:
                # keys here are used as legend labels in the outputted figure
                output_dict[str(k) + ' rep'] = float(transitions[k]) / value_sum

        if 0 in transitions:
            output_dict[str(1) + ' rep'] += float(transitions[0]) / value_sum

        datadict_array.append(output_dict)

    stacked_data = [[str(x) + ' rep' for x in replacement_order]]

    snp_sum = 0.0
    total_sum = 0.0

    for x, i in zip(aa_order, range(0, len(aa_order))):

        data = datadict_array[i]

        temp = [x]
        for y in replacement_order:

            # lumps this data into the SNP = 1 category (as it will generally be)
            if '0 rep' in data and y == 1:
                temp.append(str(data[str(y) + ' rep'] + data[str(0) + ' rep']))
            else:
                temp.append(str(data[str(y) + ' rep']))

            if y <= 1:
                snp_sum += data[str(y) + ' rep']
            total_sum += data[str(y) + ' rep']

        stacked_data.append(temp)

    glib.stackedbar([str(x) + ' rep' for x in replacement_order],
                    datadict_array, [x + '(' + str(total_aa_counter[x]) + ')' for x in aa_order],
                    'Amino Acids',
                    'AA Replacement Distance',
                    os.path.join(constants.OUTPUT_DIR, 'Estimated Nucleotide Transitions.pdf'),
                    rotation='vertical',
                    cmap='viridis')

    snp_faction = snp_sum / total_sum

    return snp_faction


def generate_heatmap(aa_order, count_dict, gene_fraction, aa_distance_dict, aa_to_codons):

    """
    
    Helper method for generating the AA => AA' transition heatmap. Count dict represents the number of AA => AA'
    occurrences that are recorded in the resistome. These values are normalized by the total number of AA changes
    in the database, and then a flexibility normalization factor that accounts for the number of possible
    codon pairs connecting two AAs and the total number of codon pairs possible.
    
    The statistical significance of the distance between AAs and their interconversion frequency is also calculated,
    and both the test statistic and p-vlaue for the Spearman correlation between the two are returned.
    
    :param aa_order: 
    :param count_dict: 
    :param gene_fraction: 
    :param aa_distance_dict: 
    :param aa_to_codons: 
    :return: 
    """

    matrix = []

    enrichment = []
    distance = []

    for key_x in aa_order:

        vector = []
        for key_y in aa_order:

            value = 0

            if key_y in count_dict[key_x]:

                # normalization factor:
                # numerator is the number of codon pairs possible (product of codon list lengths)
                # denominator: number of possible codon pairs given the genetic code
                flexibility_norm_factor = float(len(aa_to_codons[key_y]) * len(aa_to_codons[key_x])) \
                                          / (64.0 ** 2 - 64.0)

                denominator = float(sum([count_dict[key_x][v] for v in count_dict[key_x]]))

                value = float(count_dict[key_x][key_y]) / denominator * gene_fraction[key_x][
                    key_y] * 1 / flexibility_norm_factor

                vector.append(min(value, 2.0))
            else:
                vector.append(min(value, 2.0))

            enrichment.append(float(count_dict[key_x][key_y]) if key_y in count_dict[key_x] else 0)
            distance.append(aa_distance_dict[key_x][key_y])

        matrix.append(vector)

    test_stat, p_value = scipy.stats.spearmanr(enrichment, distance)

    glib.generate_heatmap(numpy.array(matrix),
                          aa_order,
                          'Normalized AA->AA\' frequency',
                          os.path.join(constants.OUTPUT_DIR,
                                       'Normalized (by unique genes) AA to AA Transition Frequency.pdf'),
                          cmap='viridis')

    return test_stat, p_value


def hamming_distance(str1, str2):

    """
    
    Calculates hamming distance between two different strings of the same length.
    
    :param str1: 
    :param str2: 
    :return: 
    """

    distance = 0
    if len(str1) != len(str2):
        raise AssertionError('unequal length in strings')

    # calculate hamming distance
    for base_x, base_y in zip(str1, str2):
        if base_x != base_y:
            distance += 1

    return distance


def codon_hamming_distance():


    """
    
    Calculates the Hamming distance between every possible codon and returns a dict (aa : dict (aa : int)) of the 
    minimum number of codon base changes required to convert AA => AA'.
    
    :return: 
    """

    code, aa_to_codon = sql_interface.get_genetic_code()

    heatmap_data = defaultdict(dict)
    output_keys = sorted(aa_to_codon.keys())
    seen = set()
    for aa_x in output_keys:
        for aa_y in output_keys:

            if aa_x == aa_y:
                heatmap_data[aa_x][aa_x] = 0
                continue

            if (aa_x, aa_y) in seen:
                continue

            seen.add((aa_x, aa_y))
            seen.add((aa_y, aa_x))

            # find the minimum of hamming distance for all pairs generated by the cartesian product of
            # codon list X and codon list Y
            cp = itertools.product(aa_to_codon[aa_x], aa_to_codon[aa_y])
            min_distance = min([hamming_distance(x[0], x[1]) for x in cp])

            heatmap_data[aa_x][aa_y] = min_distance
            heatmap_data[aa_y][aa_x] = min_distance

    return heatmap_data, aa_to_codon


def generate_aa_transition_tuples():

    """
    
    Extracts tuples of (original aa, mutated aa, position) from the Resistome.
    
    Note the different tuple ordering compared to the same method in the ProteinAnalysis file, and the lack
    of filtering to ensure the underlying sequence matches the reported WT AA.
    
    :return: 
    """

    records = sql_interface.get_gene_mutation_tuples(['aa_snps'])

    gene_count = defaultdict(dict)

    output = []

    for r in records:

        mutations = r['annotation']['aa_snps']

        for aa_mut in mutations:

            original_aa = aa_mut[1]
            modified = aa_mut[2]

            position = aa_mut[0]

            if modified not in gene_count[original_aa]:
                gene_count[original_aa][modified] = set()

            gene_count[original_aa][modified].add(r['mg1655_name'])

            output.append((original_aa, modified, position))

    return gene_count, output


def aa_chemical_property_survey(aa_tuples):

    aa_chemical_properties = dict()

    charged_aa = ['R', 'K', 'D', 'E']
    polar_aa = ['Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W']
    hydrophobic = ['A', 'I', 'L', 'M', 'F', 'V', 'P', 'G']

    for x in charged_aa:
        aa_chemical_properties[x] = 'charged'

    for x in polar_aa:
        aa_chemical_properties[x] = 'polar'

    for x in hydrophobic:
        aa_chemical_properties[x] = 'hydrophobic'

    aa_chemical_properties['*'] = 'stop'

    chemical_change_dict = defaultdict(int)

    for aa_x, aa_y, _ in aa_tuples:

        chemical_class_x = aa_chemical_properties[aa_x]
        chemical_class_y = aa_chemical_properties[aa_y]

        chemical_change_dict[(chemical_class_x, chemical_class_y)] += 1

    keys = sorted(chemical_change_dict.keys())

    for (class1, class2) in keys:

        print('Conversion between %s to %s count: %i' % (class1, class2, chemical_change_dict[(class1, class2)]))


def run_analysis():

    output_file = []

    aa_distance_dict, aa_to_codons = codon_hamming_distance()
    gene_count, output = generate_aa_transition_tuples()

    aa_chemical_property_survey(output)

    gene_fraction = compute_gene_fractions(gene_count)

    overall_aa_count, transition_count_dict, aa_transition_weighted = build_aa_transition_dict(output, aa_distance_dict)

    aa_order = sorted(transition_count_dict.keys(), reverse=True)

    temp_order = []

    for aa in aa_order:

        if aa != '*':
            temp_order.append(aa)

    temp_order.insert(0, '*')
    aa_order = temp_order

    snp_fraction = generate_aa_stacked_figure(aa_order, overall_aa_count, aa_transition_weighted)

    output_file.append('Fraction of SNP versus all transitions: %g' % (snp_fraction,))

    test_stat, p = generate_heatmap(aa_order, transition_count_dict, gene_fraction, aa_distance_dict, aa_to_codons)

    output_file.append('Spearman correlations between AA->AA transition frequency and D(i,j) between codons: %g, %g' \
          % (test_stat, p))

    with open(os.path.join(constants.OUTPUT_DIR, 'Amino Acid Change Statistical Analysis.txt'), 'w') as f:

        for line in output_file:
            f.write(line + '\n')


if __name__ == '__main__':

    run_analysis()
