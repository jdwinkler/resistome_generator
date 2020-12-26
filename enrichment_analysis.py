from resistome import constants
from resistome.sql import sql_interface
from resistome.utils import go_enrichment
import os


def gene_ontology_analysis(file_name, gene_list, enrichment_obj=None):

    if enrichment_obj is None:
        enrichment_obj = go_enrichment.EnrichmentCalculator(sql_interface.get_gene_label_relationships('go'))

    filtered_tags = enrichment_obj.compute_enrichment(gene_list,
                                                      p_value_filter=constants.P_VALUE_THRESHOLD,
                                                      apply_multih_correction=False)

    go_tag_name = sql_interface.get_go_terms([x[0] for x in filtered_tags])

    with open(file_name, 'w') as fhandle:

        for term in filtered_tags:
            fhandle.write('\t'.join((term[0], go_tag_name[term[0]], str(term[1]), str(term[2]))) + '\n')

    return filtered_tags


if __name__ == '__main__':

    genotypes = sql_interface.get_gene_mutation_tuples()
    gene_list = [record['mg1655_name'] for record in genotypes]

    print('Pulling GO terms')
    accession_go_tuples = sql_interface.get_gene_label_relationships('go')

    print('Finished getting GO terms')

    enrich_obj = go_enrichment.EnrichmentCalculator(accession_go_tuples)

    gene_ontology_analysis(os.path.join(constants.OUTPUT_DIR,
                                        'Enriched Genes (Over Database).txt'),
                           gene_list,
                           enrichment_obj=enrich_obj)

    for pcategory in sql_interface.get_phenotype_classes():

        print('Processing %s GO enrichment' % pcategory)

        genotypes = sql_interface.get_gene_mutation_tuples(filter_for_phenotype_class=[pcategory])
        gene_list = [record['mg1655_name'] for record in genotypes]

        gene_ontology_analysis(os.path.join(constants.OUTPUT_DIR,
                                            'Enriched Genes (%s phenotypes).txt' % pcategory),
                               gene_list,
                               enrichment_obj=enrich_obj)
