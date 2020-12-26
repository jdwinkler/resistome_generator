from resistome.sql import sql_interface
import os
from collections import defaultdict
import networkx
from resistome.graphics import visualization
from resistome import constants
import numpy

"""

given a mutation X, I want to know:

    1. what genes are regulated by this gene product
    2. where these genes are in the known E. coli regulatory network
    3. can we find genes that are likely to impact particular phenotypes based on their regulatory relationship?
    4. conversely, can we pin down which regulatory subset of a global regulator is responsible for a given phenotype?
        (this will also require more information from biocyc-add genome features table)

"""


def generate_matrix(gene_gene_phenotype_dict, threshold_weight=0.0):
    """
    
    Converts a dict of dicts into a matrix of ints/floats and returns the matrix, the row order, and the column
    order.
    
    :param gene_gene_phenotype_dict: 
    :param threshold_weight: 
    :return: 
    """

    cols = set()
    rows_hit = set()

    regulated_weights = defaultdict(list)
    regulated_to_regulator = defaultdict(str)

    for g1 in gene_gene_phenotype_dict:
        for g2 in gene_gene_phenotype_dict[g1]:
            weight = len(gene_gene_phenotype_dict[g1][g2])
            regulated_weights[g2].append(weight)
            regulated_to_regulator[g2] = g1

    for g2 in regulated_weights:
        if numpy.max(regulated_weights[g2]) > threshold_weight:
            cols.add(g2)
            rows_hit.add(regulated_to_regulator[g2])

    matrix = []

    rows = sorted(list(rows_hit))
    cols = sorted(list(cols))

    for g1 in rows:
        temp = []
        for g2 in cols:
            temp.append(len(gene_gene_phenotype_dict[g1].get(g2, [])))
        matrix.append(temp)

    return matrix, rows, cols


def generate_weighted_graph(gene_gene_phenotype_dict, directed=True):
    """
    
    Generated a NetworkX (un)directed graph object from a dict of dicts
    
    :param gene_gene_phenotype_dict: 
    :param directed: 
    :return: 
    """

    edges = []

    if directed:
        G = networkx.DiGraph()
    else:
        G = networkx.Graph()

    for g1 in gene_gene_phenotype_dict:
        for g2 in gene_gene_phenotype_dict[g1]:
            weight = len(gene_gene_phenotype_dict[g1][g2])
            edges.append((g1, g2, {'weight': weight}))
            if not directed:
                edges.append((g2, g1, {'weight': weight}))

    G.add_edges_from(edges)

    return G


def compare_regulatory_to_gene_mutations(mutation_data,
                                         regulatory_genes,
                                         gene_regulon_dict,
                                         directed=True,
                                         exclusions=set()):
    """
    
    Attempts to identify how regulatory mutations are linked to phenotypes as follows:
    
    1. Sorts mutation dataset into regulator/non-regulator classes (can use other binary classification schemse if you
    like); takes account directionality if directed flag = True
    
    2. Checks the intersection of the non-regulated genes leading to phenotypes with their regulators that trigger
    the same phenotypes
    
    3. Returns a dict[str: gene] : dict[str: regulator] = set of phenotype interactions
    
    :param mutation_data: 
    :param regulatory_genes: 
    :param gene_regulon_dict: 
    :param directed: 
    :param exclusions: 
    :return: 
    """

    regulator_phenotypes_specific = defaultdict(set)
    regulator_phenotypes_general = defaultdict(set)

    nonregulator_phenotypes_specific = defaultdict(set)
    nonregulator_phenotypes_general = defaultdict(set)

    for mutation in mutation_data:

        gene = mutation['mg1655_name']

        if gene in exclusions:
            continue

        if directed:
            if gene in regulatory_genes:
                regulator_phenotypes_specific[gene].add((mutation['phenotype'],
                                                         mutation['phenotype_type']))
                regulator_phenotypes_general[gene].add((mutation['phenotype_class'],
                                                        mutation['phenotype_type']))
            else:
                nonregulator_phenotypes_specific[gene].add((mutation['phenotype'],
                                                            mutation['phenotype_type']))
                nonregulator_phenotypes_general[gene].add((mutation['phenotype_class'],
                                                           mutation['phenotype_type']))
        else:
            regulator_phenotypes_specific[gene].add((mutation['phenotype'],
                                                     mutation['phenotype_type']))
            regulator_phenotypes_general[gene].add((mutation['phenotype_class'],
                                                    mutation['phenotype_type']))
            nonregulator_phenotypes_specific[gene].add((mutation['phenotype'],
                                                        mutation['phenotype_type']))
            nonregulator_phenotypes_general[gene].add((mutation['phenotype_class'],
                                                       mutation['phenotype_type']))

    potentially_causative_genes = defaultdict(dict)

    for gene in regulatory_genes:

        if gene not in gene_regulon_dict:
            continue

        regulator_phenotypes = set(regulator_phenotypes_general[gene])
        regulon = gene_regulon_dict[gene]

        for (regulated_gene, direction) in regulon:
            if regulated_gene in nonregulator_phenotypes_general:
                # which phenotypes do the mutated regulator and its regulated gene share?
                phenotype_intersection = regulator_phenotypes.intersection(set(nonregulator_phenotypes_general[regulated_gene]))
                if len(phenotype_intersection) > 0:
                    potentially_causative_genes[gene][regulated_gene] = phenotype_intersection

    return potentially_causative_genes


def generate_heatmap_plot(gg_dict, output_filename, threshold_weight=0.0):
    matrix, rows, cols = generate_matrix(gene_gene_phenotype_dict=gg_dict,
                                         threshold_weight=threshold_weight)

    visualization.generate_heatmap_asymmetric(matrix,
                                              cols,
                                              rows,
                                              'Phenotypes in Common',
                                              output_filename,
                                              x_label_rotation='vertical',
                                              dim_override=(len(cols) / 4, len(rows) / 5))


def run_analysis():
    go_regulatory_terms = ['GO:0016987', 'GO:0006351', 'GO:0006353']

    # todo: add in gene-gene network clustering

    mutation_data = sql_interface.get_gene_mutation_tuples()

    genes = set([x['mg1655_name'] for x in mutation_data])

    regulatory_genes = sql_interface.filter_genes_by_feature(genes,
                                                             go_regulatory_terms,
                                                             'go',
                                                             query_species='mg1655')

    gene_regulon_dict = sql_interface.get_gene_regulons(regulatory_genes,
                                                        query_species='mg1655')

    gene_ppi_dict = sql_interface.get_gene_physical_interactions(genes,
                                                                 query_species='mg1655')

    reg_phenotype_dict = compare_regulatory_to_gene_mutations(mutation_data,
                                                              regulatory_genes,
                                                              gene_regulon_dict)

    ppi_phenotype_dict = compare_regulatory_to_gene_mutations(mutation_data,
                                                              genes,
                                                              gene_ppi_dict,
                                                              directed=False)

    generate_heatmap_plot(reg_phenotype_dict,
                          os.path.join(constants.OUTPUT_DIR, 'GG-R-HeatMap.pdf'),
                          threshold_weight=10.0)

    generate_heatmap_plot(ppi_phenotype_dict,
                          os.path.join(constants.OUTPUT_DIR, 'GG-P-HeatMap.pdf'),
                          threshold_weight=15.0)


if __name__ == '__main__':
    run_analysis()
