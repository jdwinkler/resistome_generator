from resistome.sql import sql_interface
from resistome import constants
from resistome.utils.go_enrichment import EnrichmentCalculator
import os


def run_analysis():

    metabolite_genotypes = sql_interface.get_mutant_genotypes_as_features(type_of_feature='metabolite')

    metabolite_counts = sql_interface.get_gene_label_relationships('metabolite')

    calc = EnrichmentCalculator(metabolite_counts)

    metabolite_list = []
    for record in metabolite_genotypes:
        metabolite_list.extend(record['metabolite'])

    filtered_metabolites = calc.compute_enrichment(metabolite_list,
                                                   p_value_filter=constants.P_VALUE_THRESHOLD,
                                                   apply_multih_correction=False,
                                                   feature_conversion=False)

    with open(os.path.join(constants.OUTPUT_DIR, 'Metabolite Enrichment vs Background.txt'), 'w') as f:

        for x in filtered_metabolites:
            f.write('\t'.join(map(str, x)) + '\n')


if __name__ == '__main__':
    run_analysis()
