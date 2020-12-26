import os
from resistome import constants
from resistome.sql import sql_interface

__author__ = 'jdwinkler'


def output_circos_multigene_highlights(chromosome,
                                       mutation_data_tuples,
                                       output_file,
                                       write_file=True,
                                       append=False):

    """
    
    Generates highlight file using standard circos output format (chromosome start end annotations).
    
    Color scheme is fixed as follows:
    
    Random mutations are very very dark orange
    Overexpression mutations are red
    Deletions are blue
    Rearrangement (large_inversions) are yellow
    Indels are black
    
    Placement is as follows (from inside to outside):
    
    Random mutations are track 1
    Overexpression mutations are track 2
    Deletions are track 3
    Rearrangement (large_inversions) are track 4
    Indels are track 5
    
    :param chromosome: chromosome to place highlights on
    :param mutation_data_tuples: iterable of (int, int, string) (start, end, mutation_type) tuples
    :param output_file: output file name
    :param write_file: true/false
    :param append: append to existing file (true/false)
    :return: None
    """

    def color_scheme(mutation_type):

        if mutation_type in constants.RANDOM_mutations:
            return 'vvdorange'
        if mutation_type in constants.OE_mutations:
            return 'red'
        if mutation_type in constants.DEL_mutations:
            return 'blue'
        if mutation_type in constants.REARRANGEMENT_mutations:
            return 'yellow'
        if mutation_type in constants.INDEL_mutations:
            return 'black'

        return 'pink'

    def position_scheme(mutation_type):

        if mutation_type in constants.RANDOM_mutations:
            return 1.00, 1.06
        if mutation_type in constants.OE_mutations:
            return 1.06, 1.12
        if mutation_type in constants.DEL_mutations:
            return 1.12, 1.18
        if mutation_type in constants.REARRANGEMENT_mutations:
            return 1.18, 1.24
        if mutation_type in constants.INDEL_mutations:
            return 1.24, 1.30

        return 1.30, 1.36

    output_data = []

    for (location1, location2, mutation_type) in mutation_data_tuples:

        if location1 is None or location2 is None:
            continue

        # color = color_scheme(mutation_type)
        color = color_scheme(mutation_type)

        (r0, r1) = position_scheme(mutation_type)
        r0 = str(r0) + 'r'
        r1 = str(r1) + 'r'

        output_data.append((chromosome,
                            str(location1),
                            str(location2),
                            'fill_color=%s,stroke_color=%s,r0=%s,r1=%s' % ('d' + color,color,r0,r1)))

    if write_file:
        with open(output_file, 'w' if not append else 'a') as fhandle:
            for line in output_data:
                fhandle.write('\t'.join(line) + '\n')

    return output_data


def output_circos_linkage_annotations(network,
                                      included_genes,
                                      location_data,
                                      link_color,
                                      output_file,
                                      exclusions=set(),
                                      weight_threshold=5,
                                      weight_bonus=0,
                                      write_file=True,
                                      append=False):

    """
    
    Generate linkage information to display (gene1, gene2) and (diff exp. 1, diff exp. 2) linkages.
    
    :param network: networkx object where nodes are linked if node1, node2 is mutated or differentially expressed
    in the same host
    :param included_genes: set of genes to include 
    :param location_data: location of genes on the E. coli chromosome (dict (bnumber : (start int, end int))
    :param link_color: circos color to use for the links
    :param output_file: output file name
    :param exclusions: genes to exclude (set of str)
    :param weight_threshold: minimum occurrence to consider link for inclusion on plot
    :param weight_bonus: minimum weight of link
    :param write_file: actually write file (true/false)
    :param append: append to output_file (true/false)
    :return: 
    """

    output = []

    modified_nodes = included_genes & set(network.nodes())

    missing_set = set()

    for node in modified_nodes:

        neighbors = set(network.neighbors(node)) & included_genes

        if node in exclusions:
            continue

        if node not in missing_set:
            missing_set.add(node)

        for neighbor in neighbors:

            if neighbor in exclusions:
                continue

            if neighbor not in location_data:
                missing_set.add(neighbor)

            if network[node][neighbor]['weight'] < weight_threshold:
                continue

            weight = min(network[node][neighbor]['weight'], 5)

            # usually this is for "weird" genes that are mislabeled or are non-standard
            if node not in location_data:
                print('Gene missing from location data: %s' % (node))
                print('Skipping...')
                continue

            if neighbor not in location_data:
                print('Neighboring gene missing %s' % neighbor)
                print('Skipping')
                continue

            (chr_s, start_s, end_s) = location_data[node]
            (chr_d, start_d, end_d) = location_data[neighbor]

            output_string = chr_s + '\t' + '\t'.join([str(start_s), str(end_s)]) + '\t' + chr_d + '\t' + '\t'.join(
                [str(start_d), str(end_d)]) \
                            + '\tthickness=' + str(weight + weight_bonus) + ',color=%s' % link_color
            output.append(output_string)

    if write_file:
        with open(output_file, 'w' if not append else 'a') as fhandle:
            for line in output:
                fhandle.write(line + '\n')

    return output


def output_circos_labels(target_genes, location_data, output_file, write_file=True, append=False):

    """
    
    Outputs labels for target genes (based on number of times mutated).
    
    :param target_genes: iterable of target genes
    :param location_data: dict of (str: (int, int) where str is a gene accession, (int1 = start, int2=end)
    :param output_file: output file name
    :param write_file: actually write file (true/false)
    :param append: append to file (true/false)
    :return: 
    """

    output = []

    for tup in target_genes:
        (gene, height) = tup

        (chromosome, start, end) = location_data[gene]

        output_string = chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + gene
        output.append(output_string)

    if write_file:
        with open(output_file, 'w' if not append else 'a') as fhandle:
            for line in output:
                fhandle.write(line + '\n')

    return output


def output_circos_gene_annotations(included_genes,
                                   location_data,
                                   node_frequency,
                                   output_file,
                                   r0=1.0,
                                   label_threshold=20,
                                   write_file=True,
                                   append=False):

    """
    
    Output a "spike" to indicate how often a given gene is mutated visually. However, spikes are proportional but
    not equal to the raw mutation count for each gene to avoid the all spike no circle effect.
    
    :param included_genes: iterable of gene accessions
    :param location_data: dict of (str: (int, int) where str is a gene accession, (int1 = start, int2=end)
    :param node_frequency: number of times a given gene is mutated (dict (str: int))
    :param output_file: output file name
    :param r0: start location of the spike (float)
    :param label_threshold: minimum amount of modifications to draw spike
    :param write_file: actually write file (true/false)
    :param append: append to target file (true/false)
    :return: 
    """

    output = []

    modified_genes = set(location_data.keys()) & included_genes

    min_height = 20
    max_height = 250
    labeled_peaks = []

    for gene in modified_genes:

        (chromosome, start, end) = location_data[gene]

        height = float(node_frequency[gene]) + min_height
        height = min(height, max_height)

        # label peaks above the median threshold
        if node_frequency[gene] > label_threshold:
            labeled_peaks.append((gene, height))

        output_string = chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + gene + '\tr0=' + str(
            r0) + 'r,r1=' + str(r0) + 'r+' + str(height) + 'p'

        output.append(output_string)

    if write_file:
        with open(output_file, 'w' if not append else 'a') as fhandle:
            for line in output:
                fhandle.write(line + '\n')

    labels = output_circos_labels(labeled_peaks,
                                  location_data,
                                  output_file + '_labels.txt',
                                  write_file=write_file,
                                  append=append)

    return output, labels


def load_gene_locations(species_set, circos_chromosome_name):

    """
    
    Builds location dict for E. coli pangenome used to underlie the circos output. Location data are extracted
    from biocyc MG1655, REL606, and W databases.
    
    :param species_set: species to include (w, mg1655, rel606)
    :param circos_chromosome_name: str, name of circos chromosome (usually chr1)
    :return: dict of (str: (int, int) where str is a gene accession, (int1 = start, int2=end)
    
    """

    output_dict = dict()

    for species in species_set:

        location_data = sql_interface.get_gene_location_data(query_species=species)

        for gene in location_data:

            output_dict[gene] = (circos_chromosome_name,
                                 location_data[gene]['start'],
                                 location_data[gene]['stop'])

    return output_dict


def generate_ecoli_config(output_prefix, append=False):

    """
    
    Helper function to generate the required files for display the E. coli "mutome", with each type of mutation
    (random, deletion, amplification, inversion) displayed on a separate ring of the plot. Gene mutation frequencies
    are displayed as spikes. 
    
    :param cursor: database cursor for the resistome
    :param output_prefix: prefix to add to default output file name
    :param append: append to output files (true/false)
    :return: None
    """

    included_genes, node_frequency = sql_interface.get_gene_gene_graph(source='genetic',
                                                                       edge_count_threshold=10)

    location_data = load_gene_locations(constants.SPECIES_LIST, 'chr1')

    mutation_tuples = sql_interface.get_mutation_location_data(location_data)

    output_circos_gene_annotations(included_genes,
                                   location_data,
                                   node_frequency,
                                   os.path.join(constants.OUTPUT_DIR, output_prefix + 'ecgenes_test.txt'),
                                   r0=1.3,
                                   label_threshold=30,
                                   append=append)

    # output_circos_linkage_annotations(G,
    #                                   included_genes,
    #                                   location_data,
    #                                   'green',
    #                                   os.path.join(constants.OUTPUT_DIR, output_prefix + 'ecgenes_links.txt'),
    #                                   append=append)

    included_genes, node_frequency = sql_interface.get_gene_gene_graph(source='expression',
                                                                       edge_count_threshold=10)

    output_circos_gene_annotations(included_genes,
                                   location_data,
                                   node_frequency,
                                   os.path.join(constants.OUTPUT_DIR, output_prefix + 'ecgenes_test.txt'),
                                   r0=1.3,
                                   label_threshold=30,
                                   append=True)

    # output_circos_linkage_annotations(G,
    #                                   included_genes,
    #                                   location_data,
    #                                   'red',
    #                                   os.path.join(constants.OUTPUT_DIR, output_prefix + 'ecgenes_links.txt'),
    #                                   append=True)

    output_circos_multigene_highlights('chr1',
                                       mutation_tuples,
                                       os.path.join(constants.OUTPUT_DIR,
                                                    output_prefix + 'ecgenomic_highlights.txt'),
                                       append=append)

if __name__ == '__main__':

    generate_ecoli_config('gene_')