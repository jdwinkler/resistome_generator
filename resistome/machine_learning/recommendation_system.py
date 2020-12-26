from resistome import constants
import itertools
import numpy
import os
from resistome.graphics import visualization as glib
from resistome.sql import sql_interface
import pickle


class Feature(object):
    """
    
    Encapsulates the concept of a "feature" to enable comparison between different Resistome genotypes
    in a quasi-mathematical way.
    
    """

    def __init__(self, feature_name, feature_details):

        self.name = feature_name.upper()
        self.details = set(feature_details)

    def __eq__(self, other):

        if not isinstance(other, Feature):
            return False
        else:
            return self.name == other.name and self.details == other.details

    def __hash__(self):

        return hash(self.name)

    def __str__(self):

        return self.name


class Vector:
    """
    
    Encapsulates a set of features to enable comparison between different sets (genotypes).
    
    """

    def __init__(self, internal_id, raw_features):

        self.id = internal_id
        self.feature_set = self.build_internal_features(raw_features)

    def build_internal_features(self, feature_set):

        output_set = set()
        output_set.update(feature_set)

        return output_set

    def integer_vector(self, feature_to_int, maximum):

        """
        
        This method originally was used to encode mutational information in a numeric form, as there used
        to be a method that would convert a genotype into a binary vector (1 if a feature is mutated, 0 otherwise).
        Kept here in case there is interest.
        
        :param feature_to_int: dict (Feature: int) mapping a given feature to a list position
        :param maximum: int, length of feature vector
        :return: dict(int : int) representing the converted vector
        """

        backing = [0] * maximum

        for f in self.feature_set:

            if constants.OE_mutations >= f.details:
                backing[feature_to_int[f.name]] = 1
            elif constants.RANDOM_mutations >= f.details:
                backing[feature_to_int[f.name]] = 2
            elif constants.DEL_mutations >= f.details:
                backing[feature_to_int[f.name]] = 3
            elif constants.ADD_mutations >= f.details:
                backing[feature_to_int[f.name]] = 4
            else:
                backing[feature_to_int[f.name]] = 5

        return backing

    def __len__(self):

        return len(self.feature_set)

    def combine(self, v2):

        return Vector(self.id + v2.id, self.feature_set.union(v2.feature_set))

    def union(self, v2):

        return len(self.feature_set.union(v2.feature_set))

    def intersection(self, v2):

        return len(self.feature_set.intersection(v2.feature_set))

    def difference(self, v2):

        return len(self.feature_set.difference(v2.feature_set))

    def distance(self, v2, method='jaccard'):

        """
        
        Computes distance between two feature vectors. V2 is the vector being compared to self (V1).
        
        Similarity: number of shared eleements between V1, V2 divided by the length of V1
        Cosine similarity: value of 0 [not similar at all] to 1 [identical] based on calculating the inverse cosine
        Ochiai similarity: similar to cosine distance but normallized by the length of V1 x V2
        Jaccard similarity: amount of overlap between V1 and V2
        
        :param v2: 
        :param method: 
        :return: 
        """

        if len(self) == 0 and len(v2) == 0:
            return 1

        if method == 'similarity':
            return float(self.intersection(v2)) / float(len(self))

        if method == 'cosine':
            return 1 - float(self.intersection(v2)) / (numpy.sqrt(len(self)) * numpy.sqrt(len(v2)))

        if method == 'ochiai':
            product = set()
            for combo in itertools.product(self.feature_set, v2.feature_set):
                product.add(combo)
            return 1 - float(self.intersection(v2)) / float(len(product))

        if method == 'jaccard':
            return 1.0 - float(self.intersection(v2)) / (float(len(self) + len(v2)) - float(self.intersection(v2)))

        raise AssertionError('Requested method %s not yet implemented' % method)


def generate_vector_set(type_of_feature):
    """
    
    Converts genotypes in the resistome to feature vectors, where the features can be:
    
    Genes 'gene'
    Gene Ontology terms for the genes mutated in a given genotype 'go'
    Metabolites that are known to be perturbed upon gene deletion 'metabolite'
    
    :param type_of_feature: str, flag for feature type to use ('gene', 'go', 'metabolite')
    :return: iterable of vector objects, set of all unique feature objects discovered
    """

    if type_of_feature not in {'gene', 'go', 'metabolite'}:
        raise AssertionError('Unknown feature type: %s' % type_of_feature)

    feature_matrix = []

    all_features = set()

    genotypes = sql_interface.get_mutant_genotypes_as_features(type_of_feature)

    from collections import defaultdict

    category_dict = defaultdict(list)

    for genotype in genotypes:
        mutant_id = genotype['mutant_id']
        features = genotype[type_of_feature]

        for (phenotype_cat, ptype) in zip(genotype['pclasses'], genotype['types']):
            category_dict[mutant_id].append((phenotype_cat, ptype))

        feature_set = set([Feature(x, []) for x in features])
        all_features.update(feature_set)
        feature_matrix.append(Vector(mutant_id, feature_set))

    categories = sql_interface.get_phenotype_classes()

    return feature_matrix, all_features, categories, category_dict


def build_proposed_vector(genes, type_of_feature):
    """
    
    Converts a list of provided genes into a feature vector. Assumes that the gene names have already been
    standardized.
    
    :param genes: list of standardized gene names (should not be a set)
    :param type_of_feature: str, type of feature to convert to, ('gene', 'go', 'metabolite') as above
    :return: 
    """

    features = set()

    if type_of_feature == 'gene':

        features.update(genes)

    elif type_of_feature == 'go':

        for species in constants.SPECIES_LIST:
            features.update(sql_interface.convert_genes_to_features(genes, 'go', query_species=species))

    elif type_of_feature == 'metabolite':

        for species in constants.SPECIES_LIST:
            features.update(sql_interface.convert_genes_to_features(genes, 'metabolite', query_species=species))

    else:

        raise AssertionError('Unknown feature type: %s' % type_of_feature)

    return Vector('test', [Feature(x, []) for x in features])


def pairwise_distance_vector(proposed_vector, vector_set, method='jaccard'):
    """
    
    Computes the distance between proposed_vector and every vector in vector_set.
    
    :param proposed_vector: vector to test against vector set
    :param vector_set: iterable of vectors
    :param method: distance method to use ('jaccard', 'cosine', 'ochiai', 'similarity')
    :return: iterable of tuples (mutant.id used in resistome, distance) between proposed vector and v[mutant_id]
    """

    distances = []

    for vector in vector_set:
        distances.append((vector.id, proposed_vector.distance(vector, method=method)))

    return distances


def pairwise_distance_matrix(vector_set, method='jaccard'):
    """
    
    Computes the pairwise distance between every vector in vector_set.
    
    :param vector_set: iterable of vector objects
    :param method: str, distance method to use
    :return: dict(dict : (str: float)) representing every pairwise distance score 
    """

    seen_before = set()
    distance_distribution = []

    for v1 in vector_set:

        for v2 in vector_set:

            if (v1.id, v2.id) in seen_before or (v2.id, v1.id) in seen_before:
                continue

            seen_before.add((v1.id, v2.id))

            if len(v1) == 1 or len(v2) == 1:
                continue

            distance = v1.distance(v2, method=method)
            distance_distribution.append(distance)

    return distance_distribution


def run_analysis():
    """

    Script for running basic analysis of resistome in terms of pairwise distance between genotypes.

    :return: None
    """

    ms_matrix, _, _, _ = generate_vector_set('metabolite')

    dist_distribution = pairwise_distance_matrix(ms_matrix, method='jaccard')

    glib.histogram(dist_distribution,
                   'Distance Scores',
                   'Bin Counts',
                   os.path.join(constants.OUTPUT_DIR, 'MS Jaccard Distance Distribution.pdf'))

    go_matrix, _, _, _ = generate_vector_set('go')

    dist_distribution = pairwise_distance_matrix(go_matrix, method='jaccard')

    glib.histogram(dist_distribution,
                   'Distance Scores',
                   'Bin Counts',
                   os.path.join(constants.OUTPUT_DIR, 'GO Jaccard Distance Distribution.pdf'))

    gene_matrix, _, _, _ = generate_vector_set('gene')

    dist_distribution = pairwise_distance_matrix(gene_matrix, method='jaccard')

    glib.histogram(dist_distribution,
                   'Distance Scores',
                   'Bin Counts',
                   os.path.join(constants.OUTPUT_DIR, 'Gene Jaccard Distance Distribution.pdf'))


def generate_serialized_output():
    """
    
    Generates a serialized version of the vector sets using genes, gene ontology features. This method is primarily
    used for pre-generating these sets for use in a limited resource environment (i.e. the interface).
    
    :return: None
    """

    go_matrix, _, _, _ = generate_vector_set('go')

    gene_matrix, _, _, _ = generate_vector_set('gene')

    output_dict = {'go': go_matrix,
                   'gene': gene_matrix}

    pickle.dump(output_dict, open(os.path.join(constants.OUTPUT_DIR,
                                               'Serialized Vector Sets.obj'), 'wb'))


if __name__ == '__main__':

    run_analysis()