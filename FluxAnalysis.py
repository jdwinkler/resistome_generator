import os
import cobra
from resistome import constants
import numpy
from resistome.graphics import visualization as glib
from collections import defaultdict
from resistome.sql import sql_interface

# model obtained from: https://www.ncbi.nlm.nih.gov/pubmed/29020004

BASE_MODEL = 'iML1515.xml'


def convert_stability_data(stability_dict_tuples):
    """
    
    Builds dict (accession, position, wt_aa, mut_aa) : stability score) from INPS stability data.
    
    :param stability_dict_tuples: iterable of dicts with accession, position, wt_aa, mutant_aa, score fields where
    accession is a standardized gene id, position is an integer > 0 and < len(accession sequence), wt_aa the amino acid
    at position in the original sequence, and mut_aa the mutated AA specified in resistome records. INPS stability scores
    are floats.
    :return: dict as specified above.
    """

    output_dict = {}

    for result in stability_dict_tuples:
        accession = result['accession']
        position = result['position']
        wt_aa = result['wt_aa']
        mut_aa = result['mutant_aa']
        score = result['score']

        output_dict[(accession, position, wt_aa, mut_aa)] = score

    return output_dict


def load_cobra_model(model_file_name):
    """
    
    Wrapper function for cobra.io.read_sbml_model function. (Might want to alter models automatically after reading
    them in the future, hence the reason why this function exists).
    
    :param model_file_name: path to SBML model file.
    :return: cobra model object
    """

    return cobra.io.read_sbml_model(os.path.join(constants.INPUT_DIR,
                                                 'metabolic_models',
                                                 model_file_name))


def generate_mutated_flux_distribution(model, mutation_tuple_set):
    """
    Attempts to implement gene modifications specified for a resistome genotype (after conversion by 
    compute_flux_modifier). Once changes are implemented in the model context, the model is optimized for biomass
    formation and the objective function result, flux distribution, and gene : rxn pairings return.
    
    :param model: cobra model
    :param mutation_tuple_set: iterable of (bnumber, mutation, and flux effect-float) tuples 
    :return: genes that are modified, genes that are ignore (not found in metabolic model), float of biomass flux
    dict (str_rxn_id : float) representing flux distribution, and iterable of (gene, rxn) tuples
    """

    found_genes = 0
    missing_genes = 0

    gene_rxn_tuples = []

    with model as testbed:

        for (gene, mutation, flux_effect) in mutation_tuple_set:

            try:
                gene = testbed.genes.get_by_id(gene)
                found_genes += 1

                for reaction in gene.reactions:

                    gene_rxn_tuples.append((gene, str(reaction)))

                    if flux_effect < 1:

                        reaction.lower_bound = reaction.lower_bound * flux_effect
                        reaction.upper_bound = reaction.upper_bound * flux_effect

                    else:

                        reaction.lower_bound = reaction.lower_bound / flux_effect
                        reaction.upper_bound = reaction.upper_bound * flux_effect

            except KeyError:
                missing_genes += 1
                continue

        if found_genes > 0:

            solution = testbed.optimize()

            biomass_max_flux = solution.objective_value

            flux_dict = {key: flux_value for key, flux_value in solution.fluxes.items()}

            shadow_prices = {key: shadow_value for key, shadow_value in solution.shadow_prices.items()}

        else:

            biomass_max_flux = -1
            flux_dict = dict()
            shadow_prices = dict()

    return found_genes, missing_genes, biomass_max_flux, flux_dict, gene_rxn_tuples, shadow_prices


def compute_flux_modifier(mutation):
    """
    
    Returns a factor to multiply the flux carried by the gene product.
    
    Deletions are 0 (no flux).
    OE is doubling of flux upper limit.
    Repressing mutations halve flux.
    All other mutations are assumed to have no effect on flux (unrealistically).
    
    :param mutation: resistome mutation type
    :return: float  0 =< x <= 2.0 representing a flux modifier
    """

    if mutation in constants.INDEL_mutations or mutation in constants.DEL_mutations:

        return 0

    elif mutation in constants.OE_mutations:

        return 2

    elif mutation in constants.REP_mutations:

        return 0.5

    else:

        # these are aa changes, snps, intergenic, large indels, etc
        # unknown effects, but explored in speculative flux analysis function

        return 1


def collate_flux_estimates(model, excluded_mutation_types=set()):
    """
    Processes all genotypes in the resistome database to estimate their flux distributions as a result of their 
    mutations. The scientific validity of this is questionable; flux models are poorly constrained and will often
    yield unrealistic distributions under ideal conditions. However, the eventual idea is something along the lines of
    10.15252/msb.20167028 (Zampieri et al., 2017) in using FBA + metabolomics + shadow prices of metabolites to predict
    metabolic constraints on ALE.
    
    :param cursor: cursor to access resistome database (psycopg2)
    :param model: cobra model
    :param excluded_mutation_types: mutation types in the resistome to ignore 
    :return: dict of (rxn: list [(phenotype, phenotype class, flux_value)]) for each genotype
    """

    protein_stability_data = []

    for strain in constants.SPECIES_LIST:
        protein_stability_data.extend(sql_interface.extract_protein_stability_data('INPS',
                                                                                   query_species=strain))

    gene_position_residue_dict = convert_stability_data(protein_stability_data)

    genotypes = sql_interface.get_mutant_genotypes()

    successes = 0
    failures = 0
    counter = 0

    flux_value_dict = defaultdict(list)
    shadow_price_dict = defaultdict(list)
    affected_reactions = []

    for result in genotypes:

        gene_names = result['mg1655_names']
        annotations = result['annotations']
        mutations = result['mutations']
        phenotypes = result['phenotypes']
        pclasses = result['pclasses']

        mutation_effect_tuples = []

        for name, mutation, annotation in zip(gene_names, mutations, annotations):

            if mutation in excluded_mutation_types:
                continue

            if mutation != 'aa_snps':

                mutation_effect_tuples.append((name.lower(),
                                               mutation,
                                               compute_flux_modifier(mutation)))

            else:

                # attempts to intepret the affect of residue changes (almost certainly wrongly)

                residue_change_list = annotation['aa_snps']

                net_dG_change = 0

                for (position, wt_aa, mut_aa) in residue_change_list:
                    net_dG_change += gene_position_residue_dict.get((name, position, wt_aa, mut_aa), 0)

                # so there is no direct, simple relationship between protein stability and functionality that I know of.
                # I will just assume high energy (dG > 0) proteins have enhanced rate, the rest are reduced.

                if net_dG_change != 0:
                    modifier = net_dG_change if net_dG_change > 0 else abs(1.0 / net_dG_change)
                else:
                    modifier = 0

                mutation_effect_tuples.append((name.lower(),
                                               mutation,
                                               modifier))

        found_genes, missing_genes, biomass_objective, flux_vector, gene_rxn_tuples, shadow_prices = \
            generate_mutated_flux_distribution(model, mutation_effect_tuples)

        for flux_key in flux_vector:
            flux_value_dict[flux_key].append((phenotypes, pclasses, flux_vector[flux_key]))

        for metabolite_key in shadow_prices:
            shadow_price_dict[metabolite_key].append((phenotypes, pclasses, shadow_prices[metabolite_key]))

        counter += 1

        successes += found_genes
        failures += missing_genes

        for gene, rxn in gene_rxn_tuples:
            affected_reactions.append((phenotypes, pclasses, gene, rxn))

    return flux_value_dict, shadow_price_dict, affected_reactions, successes, failures


def overall_shadow_dict_analysis(shadow_price_dict, zero_threshold=0.1):
    """

    Computes shadow price statistics using output of collage_flux_estimates.

    :param shadow_price_dict: dict of (rxn: list [(phenotype, phenotype class, flux_value)]) for each genotype
    :param zero_threshold: value below which to consider reaction flux variance as zero
    :return: ordered list of fluxes with high variance => low variance, reaction dict of rxn_id : s/mu
    """

    prices_stat_list = []

    shadow_prices = dict()

    for metabolite in shadow_price_dict:
        # shadow prices can be both positive and negative
        vector = [x[2] for x in shadow_price_dict[metabolite]]

        mean_price = numpy.mean(vector)
        std_price = numpy.std(vector)

        prices_stat_list.append((metabolite,
                                 mean_price,
                                 std_price,
                                 std_price / mean_price if mean_price > 0 else 0))

        shadow_prices[metabolite] = std_price / mean_price if mean_price > 0 else 0

    prices_stat_list = sorted(filter(lambda x: x[2] > zero_threshold, prices_stat_list),
                              key=lambda y: y[2],
                              reverse=True)

    return prices_stat_list, shadow_prices


def category_shadow_price_dict_analysis(shadow_price_dict, zero_threshold=0.1, use_resistome_classes=True):
    """

    Repackages shadow_price_dict into category specific lists.

    :param shadow_price_dict: dict of (rxn: list [(phenotype, phenotype class, shadow_price)]) for each genotype
    :param zero_threshold: value below which to consider variance equivalent to zero
    :param use_resistome_classes: use general stress categories specified in the resistome
    :return: dict (category : sorted list of shadow_prices), dict of rxn : stdev/average abs(shadow_price)
    """

    refactored_price_dict = defaultdict(dict)

    for metabolite in shadow_price_dict:

        for (specific_class, resistome_class, price_value) in shadow_price_dict[metabolite]:

            if use_resistome_classes:
                # these classes are basically solvents_biofuels, organic acids, etc
                classes_of_interest = resistome_class
            else:
                # specific phenotypes (n-butanol resistance, etc)
                classes_of_interest = specific_class

            for class_id in classes_of_interest:

                if class_id not in refactored_price_dict[metabolite]:
                    refactored_price_dict[metabolite][class_id] = []

                refactored_price_dict[metabolite][class_id].append(abs(price_value))

    price_output_dict = defaultdict(list)
    metabolite_dict = defaultdict(dict)

    for metabolite in refactored_price_dict:

        for class_id in refactored_price_dict[metabolite]:

            mean_abs_price = numpy.mean(refactored_price_dict[metabolite][class_id])
            std_abs_price = numpy.std(refactored_price_dict[metabolite][class_id])

            price_output_dict[class_id].append((metabolite,
                                                mean_abs_price,
                                                std_abs_price,
                                                std_abs_price / mean_abs_price if mean_abs_price > 0 else 0))

            if std_abs_price < zero_threshold:
                metabolite_dict[metabolite][class_id] = 0
            else:
                metabolite_dict[metabolite][class_id] = std_abs_price / mean_abs_price

    for class_id in price_output_dict:
        prices = sorted(filter(lambda x: x[2] > zero_threshold, price_output_dict[class_id]),
                        key=lambda x: x[2],
                        reverse=True)

        price_output_dict[class_id] = prices

    return price_output_dict, metabolite_dict


def overall_flux_dict_analysis(flux_value_dict, zero_threshold=0.1):
    """
    
    Computes flux statistics using output of collage_flux_estimates.
    
    :param flux_value_dict: dict of (rxn: list [(phenotype, phenotype class, flux_value)]) for each genotype
    :param zero_threshold: value below which to consider reaction flux variance as zero
    :return: ordered list of fluxes with high variance => low variance, reaction dict of rxn_id : s/mu
    """

    fluxes = []

    reaction_dict = dict()

    for rxn in flux_value_dict:

        vector = [x[2] for x in flux_value_dict[rxn]]

        mean_flux = numpy.mean(vector)
        std_flux = numpy.std(vector)

        fluxes.append((rxn, mean_flux, std_flux, std_flux / mean_flux if mean_flux > 0 else 0))

        if std_flux < zero_threshold:
            reaction_dict[rxn] = 0
        else:
            reaction_dict[rxn] = std_flux / mean_flux

    fluxes = sorted(filter(lambda x: x[2] > zero_threshold, fluxes), key=lambda x: x[2], reverse=True)

    return fluxes, reaction_dict


def category_flux_dict_analysis(flux_value_dict, zero_threshold=0.1, use_resistome_classes=True):
    """
    
    Repackages flux_value_dict into category specific lists.
    
    :param flux_value_dict: dict of (rxn: list [(phenotype, phenotype class, flux_value)]) for each genotype
    :param zero_threshold: value below which to consider variance equivalent to zero
    :param use_resistome_classes: use general stress categories specified in the resistome
    :return: dict (category : sorted list of rxn fluxes), dict of rxn : stdev/average flux
    """

    refactored_flux_dict = defaultdict(dict)

    for rxn in flux_value_dict:

        for (specific_class, resistome_class, flux_value) in flux_value_dict[rxn]:

            if use_resistome_classes:
                # these classes are basically solvents_biofuels, organic acids, etc
                classes_of_interest = resistome_class
            else:
                # specific phenotypes (n-butanol resistance, etc)
                classes_of_interest = specific_class

            for class_id in classes_of_interest:

                if class_id not in refactored_flux_dict[rxn]:
                    refactored_flux_dict[rxn][class_id] = []

                refactored_flux_dict[rxn][class_id].append(flux_value)

    flux_output_dict = defaultdict(list)
    reaction_dict = defaultdict(dict)

    for rxn in refactored_flux_dict:

        for class_id in refactored_flux_dict[rxn]:

            mean_flux = numpy.mean(refactored_flux_dict[rxn][class_id])
            std_flux = numpy.std(refactored_flux_dict[rxn][class_id])

            flux_output_dict[class_id].append((rxn, mean_flux, std_flux, std_flux / mean_flux if mean_flux > 0 else 0))

            if std_flux < zero_threshold:
                reaction_dict[rxn][class_id] = 0
            else:
                reaction_dict[rxn][class_id] = std_flux / mean_flux

    for class_id in flux_output_dict:
        fluxes = sorted(filter(lambda x: x[2] > zero_threshold, flux_output_dict[class_id]),
                        key=lambda x: x[2],
                        reverse=True)

        flux_output_dict[class_id] = fluxes

    return flux_output_dict, reaction_dict


def summary_plot_helper_function(model, summary_stats, data_type='reaction'):
    """
    
    Converts provided summary stats into an easily plottable form; model elements are converted to their
    human-readable names encoded by the model.
    
    :param model: cobra object
    :param summary_stats: iterable of (model element object, average flux or shadow price, stdev, and stdev/average)
    tuples
    :param data_type: str flag indicating 'reaction' or 'metabolite' data 
    :return: paired iterable of element name, another iterable containing stdev)
    """

    if data_type not in {'reaction', 'metabolite'}:
        raise AssertionError('Invalid data type flag: %s' % data_type)

    output_elements = []

    for (model_element, mu, s, s_divided_by_mu) in summary_stats:

        if data_type == 'reaction':
            element_name = '%s (%s)' % (model.reactions.get_by_id(model_element).name, model_element)
        elif data_type == 'metabolite':
            element_name = model.metabolites.get_by_id(model_element).name
        else:
            raise AssertionError('Unhandled case for summary plot helper: %s' % data_type)

        output_elements.append((element_name, s))

    output_elements = sorted(output_elements, key=lambda x: x[1], reverse=True)

    return [x[0] for x in output_elements], [x[1] for x in output_elements]


def output_escher_map_data(reaction_dict, filename, append=False):
    """
    
    Writes an escher formatted reaction file for displaying fluxes.
    
    :param reaction_dict: (str: float) for rxn_ids to flux values
    :param filename: output file name
    :param append: append to filename (true/false)
    :return: None
    """

    with open(filename, 'w' if not append else 'a') as fhandle:
        for rxn in reaction_dict:
            fhandle.write(','.join([rxn, str(reaction_dict[rxn])]) + '\n')


def run_analysis():

    model = load_cobra_model(BASE_MODEL)

    flux_value_dict, shadow_price_dict, affected_rxns, _, _ = collate_flux_estimates(model)

    overall_flux_summary_stats, overall_reaction_dict = overall_flux_dict_analysis(flux_value_dict)

    overall_shadow_summary_stats, overall_shadow_dict = overall_shadow_dict_analysis(shadow_price_dict)

    output_escher_map_data(overall_reaction_dict,
                           os.path.join(constants.OUTPUT_DIR,
                                        'Escher Reaction Variation Data.csv'))

    output_escher_map_data(overall_shadow_dict,
                           os.path.join(constants.OUTPUT_DIR,
                                        'Escher Reaction Variation Data.csv'),
                           append=True)

    shadow_names, shadow_values = summary_plot_helper_function(model,
                                                               overall_shadow_summary_stats,
                                                               data_type='metabolite')

    reaction_names, reaction_values = summary_plot_helper_function(model,
                                                                   overall_flux_summary_stats,
                                                                   data_type='reaction')

    glib.bargraph(shadow_names[0:30],
                  shadow_values[0:30],
                  'Metabolite Name',
                  'Shadow price stdev/average',
                  os.path.join(constants.OUTPUT_DIR,
                               'Top Metabolite Shadow Price Variability Scores.pdf'),
                  rotation='vertical')

    glib.bargraph(reaction_names[0:30],
                  reaction_values[0:30],
                  'Reaction Name',
                  'Reaction flux stdev/average',
                  os.path.join(constants.OUTPUT_DIR,
                               'Top Reaction Flux Variability Scores.pdf'),
                  rotation='vertical')

    with open(os.path.join(constants.OUTPUT_DIR, 'Reaction Flux stdev-average table'), 'w') as f:

        for (name, value) in zip (reaction_names, reaction_values):

            name_tokens = name.split(' (')

            english_name = name_tokens[0].strip()
            cobra_name = name_tokens[1].replace('(', '').strip()

            f.write('\t'.join([cobra_name, english_name, str(value)]) + '\n')


if __name__ == '__main__':
    run_analysis()
