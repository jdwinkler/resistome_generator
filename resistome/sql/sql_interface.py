from collections import defaultdict
from resistome import constants
import psycopg2
import psycopg2.extras
import psycopg2.sql
from typing import List, Tuple

species = constants.SPECIES_LIST


class ConnectionManager:

    """
    
    Class that handles details related to connecting to the Resistome database.
    
    """

    def __init__(self, username=None, password=None, db_name='resistome'):

        self.cursor, self.connection = self.get_database_cursor(user_name=username,
                                                                password=password,
                                                                database_name=db_name)

    def get_database_cursor(self,
                            user_name,
                            password,
                            database_name):
        """
    
        Returns a connection and DictCursor object for use in querying/DB modification.
    
        Assumes a database called resistome exists, although the db credentials and name
        can be changed in inputs/db_credentials/credentials.txt.
    
        :return: 
        """

        connection = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'"
                                      % (database_name, user_name, password))
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        return cursor, connection


def get_genetic_code(query_species='mg1655', code: int = 11):
    """

    Builds two objects representing the E. coli genetic code:

    codon to AA mapping (dict (codon: AA), all str)
    AA to codon reverse mapping
    
    Query species currently does nothing, but in the event that the resistome is used to host other organisms with
    different genetic codes it will be used to return a particular species' code.

    :param query_species:
    :return: 
    """

    cursor.execute('SELECT codon, aa FROM public.genetic_code WHERE code = %s', (code,))

    code = {}
    aa_to_codon = defaultdict(list)

    for result in cursor:
        # codon is a string, not a set of codons
        code[result['codon']] = result['aa']
        aa_to_codon[result['aa']].append(result['codon'])

    return code, aa_to_codon


def get_go_terms(terms_to_extract):

    """
    
    Returns a dict of go_term (GO:XXXXXXX) : name of term.
    
    :param terms_to_extract: 
    :return: 
    """

    cursor.execute('SELECT go_term, name FROM public.go_table WHERE go_term = ANY(%s)',
                   (terms_to_extract,))

    go_tag_name = dict()
    for result in cursor:
        go_tag_name[result['go_term']] = result['name']

    return go_tag_name


def get_go_adjacency_list():
    """

    Returns an adjacency list (parent, child) edges for the GO graph encoded in the database.

    :return: 
    """

    cursor.execute('SELECT parent, child FROM public.go_adj_list')

    output = []

    for result in cursor:
        output.append((result['parent'], result['child']))

    return output


def get_gene_label_relationships(type_of_feature):

    """
    
    Builds accession : identifier list; can request accession-go, accession-metabolite, accession-uniprot relationships.
    
    :param type_of_feature: 
    :return: 
    """

    if type_of_feature not in {'go', 'metabolite', 'uniprot'}:
        raise AssertionError('Unknown type_of_feature to build gene-value dict FROM: %s' % type_of_feature)

    if type_of_feature == 'go':
        field = 'go_term'
        cursor.execute('SELECT genes.accession, go_table.go_term '
                       'FROM go_terms '
                       'INNER JOIN genes on genes.gene_id = go_terms.gene_id '
                       'INNER JOIN go_table on go_terms.go_id = go_table.go_id')
    elif type_of_feature == 'uniprot':
        field = 'region'
        cursor.execute('SELECT genes.accession, uniprot.region '
                       'FROM uniprot '
                       'INNER JOIN genes on genes.gene_id = uniprot.gene_id')
    else:
        field = 'metabolite'
        cursor.execute('SELECT genes.accession, ms_ions.ion_id as metabolite'
                       'FROM metabolomics '
                       'INNER JOIN genes on genes.gene_id = metabolomics.gene_id '
                       'INNER JOIN ms_ions on ms_ions.ion_pk = metabolomics.metabolite_id ')


    output_tuples = []

    for result in cursor:
        output_tuples.append((result['accession'], result[field]))

    return output_tuples


def get_distinct_entities(entity):

    """
    
    Returns the list of unique terms (entity = 'gene', 'go', or 'metabolite') and their count.
    
    :param entity: 
    :return: 
    """

    if entity == 'gene':
        cursor.execute('SELECT distinct on (name) name FROM resistome.mutations')
    elif entity == 'go':
        cursor.execute('SELECT distinct go_term FROM public.go_table')
    elif entity == 'metabolite':
        cursor.execute('SELECT distinct ion_id FROM public.ms_ions')
    else:
        raise AssertionError('Unknown entity type requested: %s' % entity)

    records = cursor.fetchall()

    return records, len(records)


def get_uniprot_location_data(strain, accession, protein_position):

    """
    
    
    :param strain: 
    :param accession: 
    :param protein_position: 
    :return: 
    """

    results = []

    cursor.execute('SELECT region, start, stop FROM uniprot '
                   'INNER JOIN genes on uniprot.gene_id = genes.gene_id '
                   'WHERE (genes.gene_id = ANY(SELECT gene_id from genes WHERE accession = %s) '
                   'AND start <= %s AND stop >= %s)',
                   (accession, protein_position, protein_position))

    results.extend(cursor.fetchall())

    return results


def get_aa_stability_position_value(gene, target_strain,
                                    positions_wt_mut_tuples: List[Tuple[int, str, str]],
                                    method: str,
                                    skip_errors=False):

    """
    
    :param gene: 
    :param target_strain: 
    :param positions_wt_mut_tuples
    :param method: 
    :param skip_errors: 
    :return: 
    """

    if target_strain.lower() not in species:
        raise AssertionError('Unknown species requested for protein stability: %s' % target_strain)

    output = []

    cursor.execute('select score, position, wt_aa, mutant_aa FROM protein_stability '
                   'INNER JOIN genes on protein_stability.gene_id = genes.gene_id '
                   'INNER JOIN resistome.gene_standardization on (genes.accession = resistome.gene_standardization.species_accession) '
                   'WHERE (genes.gene_id = ANY(SELECT gene_id from genes WHERE accession = %s) or resistome.gene_standardization.mg1655_accession = %s) '
                   'AND method = %s',
                   (gene, gene, method))

    pos_dict = defaultdict(set)
    pos_to_wt = dict()
    for (p, wt_aa, mut_aa) in positions_wt_mut_tuples:
        pos_dict[p].add(mut_aa)
        pos_to_wt[p] = wt_aa

    output = dict()
    for x in cursor:
        if x['position'] in pos_to_wt and x['mutant_aa'] in pos_dict[x['position']]:
            x_position = x['position']
            if skip_errors and x['wt_aa'] != pos_to_wt[x_position]:
                continue

            output[(x_position, pos_to_wt[x_position], x['mutant_aa'])] = x['score']

    return output

    #
    # if skip_errors:
    #     # less stringent, will succeed unless gene is missing entirely
    #     cursor.execute('SELECT score, position FROM protein_stability '
    #                    'INNER JOIN genes on protein_stability.gene_id = genes.gene_id '
    #                    'INNER JOIN resistome.gene_standardization on (genes.accession = resistome.gene_standardization.species_accession) '
    #                    'WHERE (genes.gene_id = ANY(SELECT gene_id from genes WHERE accession = %s) or resistome.gene_standardization.mg1655_accession = %s) '
    #                    'AND position = %s '
    #                    'AND mutant_aa = %s '
    #                    'AND method = %s',
    #                    (gene, gene, position, mut_aa, method))
    # else:
    #     # forces query to only match if WT AA = that in sequence record
    #     cursor.execute('SELECT score FROM protein_stability '
    #                    'INNER JOIN genes on protein_stability.gene_id = genes.gene_id '
    #                    'INNER JOIN resistome.gene_standardization on (genes.accession = resistome.gene_standardization.species_accession) '
    #                    'WHERE (genes.accession = ANY(SELECT gene_id from genes WHERE accession = %s) or resistome.gene_standardization.mg1655_accession = %s) '
    #                    'AND position = %s '
    #                    'AND wt_aa = %s '
    #                    'AND mutant_aa = %s '
    #                    'AND method = %s',
    #                    (gene, gene, position, wt_aa, mut_aa, method))
    #
    # results = cursor.fetchone()


def get_nearby_genes_FROM_position(position, query_species='mg1655'):

    if query_species not in species:
        raise AssertionError('Unknown species requested: %s' % species)

    # get nearest gene before and after position

    cursor.execute('SELECT accession FROM genes '
                   'INNER JOIN gene_locations on gene_locations.gene_id = genes.gene_id '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'WHERE stop < %s and strain.strain = %s'
                   'order by stop DESC LIMIT 1', (position, query_species))

    nearest_gene_before = cursor.fetchone()['accession']

    cursor.execute('SELECT accession FROM genes '
                   'INNER JOIN gene_locations on gene_locations.gene_id = genes.gene_id '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'WHERE start > %s and strain.strain = %s'
                   'order by start DESC LIMIT 1', (position, query_species))

    nearest_gene_after = cursor.fetchone()['accession']

    return nearest_gene_before, nearest_gene_after


def extract_protein_stability_data(method, query_species='mg1655'):
    """

    Extracts all gene stability data for all species and returns it as a list of dicts.

    accession
    position
    wt_aa
    mut_aa
    score

    :param method: 
    :param query_species
    :return: 
    """

    if method not in {'SNAP2', 'INPS'}:
        raise AssertionError('Unknown method requested: %s' % method)

    if query_species not in species:
        raise AssertionError('Unknown species requested: %s' % species)

    final_query = 'SELECT * FROM protein_stability ' \
                  'INNER JOIN genes on genes.gene_id = protein_stability.gene_id ' \
                  'INNER JOIN strain on strain.strain_id = genes.strain_id ' \
                  'WHERE method = upper(%s) and strain.strain = %s'

    cursor.execute(final_query, (method, query_species))

    results = cursor.fetchall()

    if results is None:
        raise ValueError('No results obtained')
    else:
        return results


def get_gene_identifiers(query_species='mg1655'):

    """
    
    Returns distinct gene names for the species requested.
    
    :param query_species: 
    :return: 
    """

    if query_species not in species:
        raise AssertionError('Unknown species gene set requested: %s' % query_species)

    cursor.execute('SELECT distinct UPPER(accession) as accession FROM genes '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'WHERE strain.strain = %s', (query_species,))

    results = cursor.fetchall()

    return [x['accession'] for x in results]


def get_gene_common_names(gene_set, query_species='mg1655'):
    """

    Returns the commonly used name for the requested genes.

    :param query_species: 
    :return: 
    """

    if query_species not in species:
        raise AssertionError('Unknown species gene set requested: %s' % query_species)

    cursor.execute('SELECT upper(accession) as accession, name '
                   'FROM genes '
                   'INNER JOIN strain on genes.strain_id = strain.strain_id '
                   'WHERE accession = ANY(%s) and strain.strain = %s',
                   (list(gene_set), query_species))

    results = cursor.fetchall()

    output_dict = dict()

    for result in results:
        output_dict[result['accession']] = result['name']

    return output_dict


def get_genome_sequence(query_species='mg1655'):
    """

    Returns the genome sequence for the requested species. Does not handle multiple chromosomes.

    :param query_species: 
    :return: 
    """

    if query_species not in species:
        raise AssertionError('Unknown species genome requested: %s' % query_species)

    cursor.execute('SELECT sequence FROM genome '
                   'INNER JOIN strain on strain.strain_id = genome.strain_id '
                   'WHERE strain.strain = %s', (query_species,))

    return cursor.fetchone()['sequence']


def get_gene_location_data(query_species='mg1655'):
    """

    Returns a dict (accession : dict (start, stop: int) representing the position of all genes in the requested genome.

    :param query_species:
    :return: 
    """

    genome_data = {}

    if query_species not in species:
        raise AssertionError('Unknown species gene data requested: %s' % query_species)

    cursor.execute('SELECT start, stop, accession FROM genes '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'INNER JOIN gene_locations on genes.gene_id = gene_locations.gene_id '
                   'WHERE strain.strain = %s', (query_species,))

    for result in cursor.fetchall():
        accession = result['accession']
        if accession in genome_data:
            genome_data[accession] = {'start': min(genome_data[accession]['start'], result['start']),
                                      'stop': max(genome_data[accession]['stop'], result['stop'])}
        else:
            genome_data[accession] = {'start': result['start'],
                                      'stop': result['stop']}

    return genome_data


def get_polypeptide_sequence(accession, query_species='mg1655'):

    """
    
    Returns the polypeptide sequence associated with accession. None is returned for sequences without translated
    counterparts.
    
    :param accession: 
    :param query_species: 
    :return: 
    """

    if query_species not in species:
        raise AssertionError('Unknown species coding sequence requested: %s' % query_species)

    cursor.execute('SELECT aa_seq FROM aa_sequences '
                   'INNER JOIN genes on aa_sequences.gene_id = genes.gene_id '
                   'WHERE genes.accession = %s', (accession,))

    # don't assume any particular form of the cursor (dict, normal, etc...)
    # let the caller handle
    # quits once any result is found, since gene ids are unique
    results = cursor.fetchone()

    if results is None:
        return None
    else:
        return results['aa_seq']


def filter_genes_by_feature(genes, features, type_of_feature, query_species='mg1655'):

    features = list(features)
    genes = list(genes)

    if type_of_feature == 'go':
        cursor.execute('SELECT accession FROM go_terms '
                       'INNER JOIN genes on genes.gene_id = go_terms.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'WHERE go_term = ANY(%s) '
                       'AND gene.accession = ANY(%s) '
                       'AND strain.strain = %s', (features, genes, query_species))
    elif type_of_feature == 'metabolite':
        cursor.execute('SELECT accession FROM metabolomics '
                       'INNER JOIN genes on genes.gene_id = go_terms.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'INNER JOIN ms_ions on ms_ions.ion_pk = metabolomics_id '
                       'WHERE ms_ions.name = ANY(%s) '
                       'AND accession = ANY(%s) '
                       'AND strain.strain = %s', (features, genes, query_species))
    elif type_of_feature == 'gene':
        raise AssertionError('Cannot filter gene list by gene features!')
    else:
        raise AssertionError('Unknown feature type: %s' % features)

    return set([x['accession'] for x in cursor])


def get_gene_physical_interactions(genes, query_species='mg1655'):

    return get_gene_interactions(genes, 'ecolinet_v1', query_species=query_species)


def get_gene_regulons(genes, query_species='mg1655'):

    return get_gene_interactions(genes, 'regulondb', query_species=query_species)


def get_gene_interactions(genes, method, query_species='mg1655'):

    genes = list(genes)

    cursor.execute('SELECT '
                   'regulator as regulator, '
                   'array_agg(target) as regulon, '
                   'array_agg(direction) as regulation_direction '
                   'FROM interactions '
                   'INNER JOIN strain on strain.strain_id = interactions.strain_id '
                   'WHERE regulator = ANY(%s) and db_source = %s '
                   'AND strain.strain = %s '
                   'group by regulator', (genes, method, query_species))

    output_dict = dict()
    for result in cursor:
        output_dict[result['regulator']] = [(gene, dir) for gene, dir in zip(result['regulon'],
                                                                             result['regulation_direction'])]

    return output_dict


def get_mutant_genotypes():

    """
    
    Aggregate resistome data back into genotypes.
    
    :return: 
    """

    cursor.execute('SELECT '
                   'resistome.mutants.mutant_id, '
                   'resistome.mutants.paper_id, '
                   'resistome.mutants.strain,'
                   'array_agg(resistome.phenotypes.phenotype) as phenotypes,'
                   'array_agg(resistome.phenotypes.phenotype_class) as pclasses,'
                   'array_agg(resistome.phenotypes.phenotype_type) as types '
                   'FROM resistome.mutants '
                   'left join resistome.phenotypes on (resistome.phenotypes.mutant_id = resistome.mutants.mutant_id) '
                   'group by resistome.mutants.mutant_id')

    phenotype_results = cursor.fetchall()

    cursor.execute('SELECT '
                   'resistome.mutations.mutant_id, '
                   'array_agg(resistome.mutations.gene_id) as gene_ids, '
                   'array_agg(resistome.mutations.name) as species_names, '
                   'array_agg(resistome.gene_standardization.mg1655_accession) as mg1655_names, '
                   'array_agg(resistome.annotations.mutation_type) as mutations, '
                   'array_agg(resistome.annotations.annotation) as annotations '
                   'FROM resistome.mutations '
                   'left join resistome.annotations on (resistome.annotations.gene_id = resistome.mutations.gene_id) '
                   'INNER JOIN resistome.gene_standardization on (resistome.gene_standardization.species_accession = resistome.mutations.name '
                   '                                              AND upper(resistome.gene_standardization.strain) = upper(resistome.mutations.strain)) '
                   'group by resistome.mutations.mutant_id')

    genotype_results = cursor.fetchall()

    genotype_dicts = defaultdict(dict)

    for result in genotype_results:

        genotype_dicts[result['mutant_id']]['mutant_id'] = result['mutant_id']
        genotype_dicts[result['mutant_id']]['gene_ids'] = result['gene_ids']
        genotype_dicts[result['mutant_id']]['species_names'] = set(result['species_names'])
        genotype_dicts[result['mutant_id']]['mg1655_names'] = set(result['mg1655_names'])
        genotype_dicts[result['mutant_id']]['paired_species_names'] = result['species_names']
        genotype_dicts[result['mutant_id']]['paired_mg1655_names'] = result['mg1655_names']

        if None in result['mutations']:
            genotype_dicts[result['mutant_id']]['mutations'] = []
        else:
            genotype_dicts[result['mutant_id']]['mutations'] = result['mutations']

        if None in result['annotations']:
            genotype_dicts[result['mutant_id']]['annotations'] = []
        else:
            genotype_dicts[result['mutant_id']]['annotations'] = result['annotations']

    for result in phenotype_results:

        if 'mutant_id' not in genotype_dicts[result['mutant_id']]:
            genotype_dicts[result['mutant_id']]['mutant_id'] = result['mutant_id']

            genotype_dicts[result['mutant_id']]['species_names'] = []
            genotype_dicts[result['mutant_id']]['mg1655_names'] = []
            genotype_dicts[result['mutant_id']]['paired_species_names'] = []
            genotype_dicts[result['mutant_id']]['paired_mg1655_names'] = []
            genotype_dicts[result['mutant_id']]['mutations'] = []
            genotype_dicts[result['mutant_id']]['annotations'] = []
            genotype_dicts[result['mutant_id']]['gene_ids'] = []

        genotype_dicts[result['mutant_id']]['paper_id'] = result['paper_id']
        genotype_dicts[result['mutant_id']]['strain'] = result['strain']
        genotype_dicts[result['mutant_id']]['phenotypes'] = result['phenotypes']
        genotype_dicts[result['mutant_id']]['pclasses'] = result['pclasses']
        genotype_dicts[result['mutant_id']]['types'] = result['types']

    return genotype_dicts.values()


def get_gene_gene_graph(source='genetic', edge_count_threshold=50):

    """

    Build interaction graph to connect genes that are either mutated in the same strain or differentially
    expressed in the same strain.

    :param source: str, either 'genetic' or 'expression' to indicate what type of interaction graph is deisred
    :param edge_count_threshold: int, minimum count of edges to include
    :return: networkx interaction graph G, set of nodes included, dict (str: int) of node counts
    """

    if source == 'genetic':

        cursor.execute('SELECT '
                       'array_agg(resistome.gene_standardization.mg1655_accession) as names,'
                       'array_agg(resistome.annotations.mutation_type) as mtypes '
                       'FROM resistome.mutations '
                       'INNER JOIN resistome.annotations on '
                       '(resistome.mutations.gene_id = resistome.annotations.gene_id) '
                       'INNER JOIN resistome.gene_standardization on '
                       '(resistome.mutations.name = resistome.gene_standardization.species_accession) '
                       'GROUP BY resistome.mutations.mutant_id')

    elif source == 'expression':

        cursor.execute('SELECT '
                       'array_agg(resistome.gene_standardization.mg1655_accession) as names,'
                       'array_agg(resistome.expressions.fold_change) as mtypes '
                       'FROM resistome.expressions '
                       'INNER JOIN resistome.gene_standardization on '
                       '(resistome.expressions.name = resistome.gene_standardization.species_accession) '
                       'GROUP BY resistome.expressions.mutant_id')

    else:
        raise AssertionError('Unknown source data for graph generation')

    node_mutations = defaultdict(set)
    node_counter = defaultdict(int)
    edge_weights = defaultdict(int)
    edge_set = set()
    node_set = set()

    design_count = 0

    for result in cursor:

        names = result['names']
        m_types = result['mtypes']

        nodes_to_process = []

        for name, mutation in zip(names, m_types):

            if source == 'genetic':

                node_counter[name] += 1
                node_mutations[name].add(mutation)
                node_set.add(name)
                nodes_to_process.append(name)

            elif source == 'expression':

                # filter for weakly differentially expressed genes
                if mutation is None or abs(mutation) < 5.0:
                    continue

                node_counter[name] += 1
                node_mutations[name].add(mutation)
                node_set.add(name)
                nodes_to_process.append(name)

    return node_set, node_counter


def get_mutant_genotypes_as_features(type_of_feature):

    if type_of_feature not in {'gene', 'go', 'metabolite'}:
        raise AssertionError('Unknown feature type requested: %s' % type_of_feature)

    records = get_mutant_genotypes()
    output_results = []
    for result in records:

        genes = result['mg1655_names']

        if type_of_feature == 'go' or type_of_feature == 'metabolite':
            converted_features = convert_genes_to_features(list(genes),
                                                           type_of_feature,
                                                           query_species=result['strain'])
        else:
            converted_features = genes

        result[type_of_feature] = set(converted_features)
        output_results.append(result)

    return output_results


def get_mutation_location_data(location_data, debug_flag=False):

    """
    Extracts mutational data FROM the resistome and converts location data to absolute positioning for use by
    circos.

    :param location_data: dict of (str: (int, int) WHERE str is a gene accession, (int1 = start, int2=end)
    :param debug_flag
    :return: iterable of tuples (start, end, resistome mutation type)
    """

    cursor.execute('SELECT resistome.mutations.name, '
                   'array_agg(resistome.annotations.annotation) as annotations,'
                   'array_agg(resistome.annotations.mutation_type) as mutation_types '
                   'FROM resistome.mutations '
                   'INNER JOIN resistome.annotations on (resistome.mutations.gene_id = resistome.annotations.gene_id) '
                   'group by resistome.mutations.gene_id')

    if debug_flag:
        return cursor.fetchall()

    output_tuples = []

    for result in cursor:

        name = result['name']
        mutation_types = result['mutation_types']
        annotations = result['annotations']

        location_info = location_data.get(name, None)

        # some mutations are not specifically located in the target papers
        if location_info is None:
            continue

        if len(location_data) > 2:
            start = location_info[1]
            stop = location_info[2]
        else:
            start = location_info[0]
            stop = location_info[1]

        for mutation_type, annotation in zip(mutation_types, annotations):

            if mutation_type in constants.DEL_mutations:

                output_tuples.append((start,
                                      stop,
                                      mutation_type))

            elif mutation_type in constants.OE_mutations:

                output_tuples.append((start, stop, mutation_type))

            elif mutation_type in constants.INDEL_mutations:

                # this is skipped because IS insertions are always associated with
                # indels giving the location/size of the
                # mutation.

                if mutation_type == 'is_insertion' and 'indel' not in annotation:
                    output_tuples.append((start, stop, mutation_type))
                else:
                    for (position, size, indel_type) in annotation[mutation_type]:
                        # the locations of some mutations are textually described
                        # skip those
                        if type(position) != int or type(size) != int:
                            continue
                        output_tuples.append((position,
                                              position + size,
                                              mutation_type))

            # aa_snps, nuc_snps, etc
            elif mutation_type in constants.RANDOM_mutations:

                if mutation_type == 'nuc_snps':
                    for (location, _, _) in annotation[mutation_type]:
                        if location is None:
                            continue
                        output_tuples.append((location, location + 1, mutation_type))
                elif mutation_type == 'aa_snps':
                    for (aa_location, _, _) in annotation[mutation_type]:
                        if aa_location is None:
                            continue
                        output_tuples.append((start + aa_location * 3, start + aa_location * 3 + 3, mutation_type))

            # large_deletions, large_inversions, large_amplifications
            elif mutation_type in constants.REARRANGEMENT_mutations:

                output_tuples.append((start, stop, mutation_type))

    return output_tuples


def get_gene_mutation_tuples(filter_for_mutation_types=None,
                             filter_for_phenotype_class=None):

    """
    
    :param filter_for_mutation_types: 
    :param filter_for_phenotype_class: 
    :return: 
    """

    genotypes = get_mutant_genotypes()

    # denormalize into gene-mutation-phenotype_class

    output_dicts = []

    for genotype in genotypes:

        gene_ids = genotype['gene_ids']
        paired_species = genotype['paired_species_names']
        paired_mg1655 = genotype['paired_mg1655_names']

        mutation_types = genotype['mutations']
        annotations = genotype['annotations']
        pclasses = genotype['pclasses']

        if filter_for_phenotype_class is not None and set(pclasses).isdisjoint(filter_for_phenotype_class):
            continue

        for (gene_id, species_name, mg1655_name, mutation, annotation) in zip(gene_ids, paired_species,
                                                                              paired_mg1655, mutation_types,
                                                                              annotations):

            if filter_for_mutation_types is not None and mutation not in filter_for_mutation_types:
                continue

            output_dicts.append({'species_name': species_name,
                                 'mg1655_name': mg1655_name,
                                 'mutation_type': mutation,
                                 'annotation': annotation,
                                 'gene_id': gene_id,
                                 'strain': genotype['strain']})

    return output_dicts


def get_paper_doe_types():

    cursor.execute('SELECT internal_name as method FROM resistome.term_explanation WHERE term_type = %s', ('DOE',))

    output_list = []
    for x in cursor:
        output_list.append(x['method'])
    return output_list


def get_mutation_counts_by_doe_type(types_to_include):

    types_to_include = list(types_to_include)

    cursor.execute('SELECT '
                   'count(resistome.annotations.mutation_type) as count '
                   'FROM resistome.mutations '
                   'INNER JOIN resistome.annotations on (resistome.annotations.gene_id = resistome.mutations.gene_id) '
                   'INNER JOIN resistome.papers on (resistome.papers.paper_id = resistome.mutations.paper_id) '
                   'WHERE resistome.papers.methods <@ %s '
                   'group by (resistome.mutations.mutant_id)',
                   (types_to_include,))

    mutation_type_counts = cursor.fetchall()

    cursor.execute('SELECT '
                   'count(distinct resistome.mutations.name) as count '
                   'FROM resistome.mutations '
                   'INNER JOIN resistome.papers on (resistome.papers.paper_id = resistome.mutations.paper_id) '
                   'WHERE resistome.papers.methods <@ %s '
                   'group by (resistome.mutations.mutant_id)',
                   (types_to_include,))

    gene_counts = cursor.fetchall()

    return [x['count'] for x in mutation_type_counts], [x['count'] for x in gene_counts]


def mutation_go_tags_by_doe(types_to_include):

    types_to_include = list(types_to_include)

    query = 'SELECT count(distinct go_terms.go_id) ' \
            'FROM resistome.mutations ' \
            'INNER JOIN genes on (resistome.mutations.name = genes.accession) ' \
            'INNER JOIN go_terms on go_terms.gene_id = genes.gene_id ' \
            'INNER JOIN strain on strain.strain_id = genes.strain_id ' \
            'INNER JOIN resistome.papers on (resistome.papers.paper_id = resistome.mutations.paper_id) ' \
            'WHERE resistome.papers.methods <@ %s AND strain.strain = ANY(%s) ' \
            'group by resistome.mutations.mutant_id '

    cursor.execute(query, (types_to_include, list(species)))

    return [x['count'] for x in cursor]


def gene_essentiality_by_doe(types_to_include):

    types_to_include = list(types_to_include)

    query = 'SELECT SUM(CASE WHEN genes.essential = True THEN 1 ELSE 0 END) as recorded_essential ' \
            'FROM resistome.mutations ' \
            'INNER JOIN genes on (resistome.mutations.name = genes.accession) ' \
            'INNER JOIN strain on strain.strain_id = genes.strain_id ' \
            'INNER JOIN resistome.papers on (resistome.papers.paper_id = resistome.mutations.paper_id) ' \
            'WHERE resistome.papers.methods <@ %s AND strain.strain = %s ' \
            'group by resistome.mutations.mutant_id '

    cursor.execute(query, (types_to_include, 'mg1655'))

    return [x['recorded_essential'] for x in cursor]


def convert_features_to_genes(features, type_of_feature, query_species='mg1655'):

    if query_species not in species:
        raise AssertionError('Unknown species: %s')

    if type_of_feature == 'go':
        cursor.execute('SELECT genes.accession FROM go_terms '
                       'INNER JOIN genes on go_terms.gene_id = genes.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'INNER JOIN go_table on go_table.go_id = go_terms.go_id '
                       'WHERE go_table.go_term = ANY(%s) and strain.strain = %s', (features, query_species))
    elif type_of_feature == 'metabolite':
        cursor.execute('SELECT genes.accession FROM metabolomics '
                       'INNER JOIN genes on metabolomics.gene_id = genes.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'INNER JOIN ms_ions on ms_ions.ion_pk = metabolomics.metabolite_id '
                       'WHERE ms_ions.name = ANY(%s) AND strain.strain = %s', (features, query_species))
    elif type_of_feature == 'gene':
        return features
    else:
        raise AssertionError('Unknown feature type: %s' % features)

    return set([x['accession'] for x in cursor])


def convert_genes_to_features(genes, type_of_feature, query_species='mg1655'):

    if query_species not in species:
        raise AssertionError('Unknown species: %s')

    if type_of_feature == 'go':
        cursor.execute('SELECT go_table.go_term as feature FROM go_terms '
                       'INNER JOIN genes on genes.gene_id = go_terms.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'INNER JOIN go_table on go_table.go_id = go_terms.go_id '
                       'WHERE genes.accession = ANY(%s) and strain.strain = %s', (genes, query_species))
    elif type_of_feature == 'metabolite':
        cursor.execute('SELECT ms_ions.ion_id as feature FROM metabolomics '
                       'INNER JOIN ms_ions on ms_ions.ion_pk = metabolomics.metabolite_id '
                       'INNER JOIN genes on genes.gene_id = metabolomics.gene_id '
                       'INNER JOIN strain on strain.strain_id = genes.strain_id '
                       'WHERE genes.accession = ANY(%s) AND strain.strain = %s', (genes, query_species))
    elif type_of_feature == 'gene':
        return genes
    else:
        raise AssertionError('Unknown feature type: %s' % genes)

    return set([x['feature'] for x in cursor])


def convert_phenotype_to_class(phenotype_name):
    """

    Converts a specific phenotype to a general phenotype class (solvents_biofuels, etc).

    :param phenotype_name: 
    :return: 
    """

    cursor.execute('SELECT specific_classes FROM resistome.phenotype_standardization '
                   'WHERE (resistome.phenotype_standardization.standard_name = %s)',
                   (phenotype_name,))

    return cursor.fetchone()


def get_phenotype_classes():

    """
    
    :return: 
    """

    cursor.execute('SELECT internal_name FROM resistome.term_explanation '
                   'WHERE upper(term_type) = upper(%s)', ('Tags',))

    categories = set([x['internal_name'] for x in cursor])

    return categories


def convert_full_identifiers_to_abbreviations(abbrev_type, keys):
    """

    Converts a list of keys to standardized representations labeled as abbrev_type.

    Allowed types: 

    categories
    journal
    methods
    mutations
    phenotype

    :param abbrev_type: 
    :param keys: 
    :return: 
    """

    filtered_keys = []
    for o_name in keys:

        cursor.execute(
            'SELECT converted_entry FROM resistome.abbreviations '
            'WHERE upper(entry_type) = upper(%s) and upper(entry) = UPPER(%s)',
            (abbrev_type, o_name))

        c_name = cursor.fetchone()

        if c_name is None:
            filtered_keys.append(o_name)
        else:
            filtered_keys.append(c_name['converted_entry'])

    return filtered_keys


def distribution_sql(schema, table, field, threshold):

    cursor.execute('SELECT %s FROM %s.%s' % (field, schema, table))
    count_dict = defaultdict(int)
    value_array = []

    for result in cursor:
        count_dict[result[field].upper()] += 1
        value_array.append(result[field])

    real_journal_dict = defaultdict(int)

    other_value = 0

    for key in count_dict:
        if count_dict[key] <= threshold:
            other_value += count_dict[key]
        else:
            real_journal_dict[key] += count_dict[key]

    key_vals = [(x, y) for x, y in real_journal_dict.items()]
    key_vals = sorted(key_vals, key=lambda x: x[1], reverse=True)
    keys = [x[0] for x in key_vals]
    keys.append('other')
    values = [y[1] for y in key_vals]
    values.append(other_value)

    return keys, values, value_array


def get_table_statistics():

    year_papers = defaultdict(int)
    year_mutation = defaultdict(list)
    year_models = defaultdict(list)

    method_counts = defaultdict(int)
    mutation_counts = defaultdict(int)

    mutant_count_array = []
    mutation_count_array = []
    total_models = []

    gene_usage = defaultdict(int)

    tolerance_phenotypes = defaultdict(int)

    specific_phenotypes = set()

    cursor.execute('SELECT resistome.papers.paper_id, '
                   'array_agg(resistome.paper_tags.tag) as tags,'
                   'resistome.papers.year, '
                   'resistome.papers.designs, '
                   'count(resistome.mutants.mutant_id) as mutant_count,'
                   'array_agg(resistome.mutants.mutant_id) as mutant_ids '
                   'FROM resistome.papers '
                   'INNER JOIN resistome.mutants on (resistome.papers.paper_id = resistome.mutants.paper_id) '
                   'INNER JOIN resistome.paper_tags on (resistome.papers.paper_id = resistome.paper_tags.paper_id) '
                   'group by resistome.papers.paper_id')

    for record in cursor.fetchall():

        for tag in record['tags']:
            tolerance_phenotypes[tag] += 1

        year_papers[record['year']] += 1
        total_models.append(record['designs'])
        year_models[record['year']].append(record['mutant_count'])
        mutant_count_array.append(record['mutant_count'])

        cursor.execute('SELECT distinct on (mutant_id, gene_id) '
                       'gene_id, '
                       'name '
                       'FROM resistome.mutations '
                       'WHERE mutant_id = ANY(%s)', (record['mutant_ids'],))

        gene_results = cursor.fetchall()
        mutation_count_array.append(len(gene_results))
        year_mutation[record['year']].append(len(gene_results))

        for gene in gene_results:
            gene_usage[gene['name']] += 1

    paper_count = sum(year_papers.values())

    cursor.execute('SELECT phenotype FROM resistome.phenotypes')
    specific_phenotypes.update([x['phenotype'] for x in cursor])

    cursor.execute('SELECT method FROM resistome.mutant_methods')
    for method in [x['method'] for x in cursor]:
        method_counts[method] += 1

    cursor.execute('SELECT count(*) as model_count FROM resistome.mutants')

    num_ecoli_designs = cursor.fetchone()['model_count']

    cursor.execute('SELECT distinct on (mutant_id, gene_id) name '
                   'FROM resistome.mutations '
                   'WHERE upper(species) = upper(%s)',
                   ('Escherichia coli',))
    num_ecoli_genes = len(cursor.fetchall())

    cursor.execute('SELECT mutation_type FROM resistome.annotations')
    results = cursor.fetchall()
    for info in results:
        mutation_counts[info['mutation_type']] += 1

    cursor.execute('SELECT mutant_id, count(*) FROM resistome.mutations '
                   'group by mutant_id')
    design_results = cursor.fetchall()

    return (year_papers,
            year_models,
            year_mutation,
            paper_count,
            specific_phenotypes,
            method_counts,
            num_ecoli_designs,
            num_ecoli_genes,
            design_results,
            gene_usage,
            tolerance_phenotypes,
            mutant_count_array,
            mutation_count_array,
            mutation_counts)


def get_aggregate_table_data():

    cursor.execute('SELECT '
                   'array_agg(distinct resistome.mutations.name) as species_name, '
                   'array_agg(distinct resistome.gene_standardization.mg1655_accession) as mg1655_name, '
                   'array_agg(distinct resistome.annotations.mutation_type) as mutation_type, '
                   'array_agg(distinct resistome.annotations.annotation) as annotation, '
                   'array_agg(resistome.phenotypes.phenotype) as phenotype, '
                   'array_agg(resistome.phenotypes.phenotype_type) as phenotype_type '
                   'FROM resistome.mutations '
                   'INNER JOIN resistome.annotations on (resistome.mutations.gene_id=resistome.annotations.gene_id) '
                   'INNER JOIN resistome.phenotypes on (resistome.mutations.mutant_id=resistome.phenotypes.mutant_id) '
                   'INNER JOIN resistome.gene_standardization on '
                   '(resistome.gene_standardization.species_accession = resistome.mutations.name) '
                   'group by resistome.mutations.gene_id')

    return cursor.fetchall()


def get_tag_statistics():

    category_counter = defaultdict(dict)
    frequency = defaultdict(int)

    cursor.execute('SELECT resistome.papers.year,'
                   'resistome.paper_tags.tag '
                   'FROM resistome.papers '
                   'INNER JOIN resistome.paper_tags on (resistome.paper_tags.paper_id = resistome.papers.paper_id)')

    for result in cursor:

        frequency[result['tag']] += 1
        if result['tag'] in category_counter[result['year']]:
            category_counter[result['year']][result['tag']] += 1
        else:
            category_counter[result['year']][result['tag']] = 1

    return frequency, category_counter


def journal_method_design_distributions(
                                        journal_threshold=5,
                                        method_threshold=10,
                                        mutation_type_threshold=20):
    journal_keys, journal_values, _ = distribution_sql(
                                                       'resistome',
                                                       'papers',
                                                       'journal',
                                                       journal_threshold)

    journal_keys = convert_full_identifiers_to_abbreviations('journal', journal_keys)

    # methods used for mutant design figure
    method_keys, method_values, _ = distribution_sql(
                                                     'resistome',
                                                     'mutant_methods',
                                                     'method',
                                                     method_threshold)
    method_keys = convert_full_identifiers_to_abbreviations('methods', method_keys)

    # mutations made figure
    mutation_keys, mutation_values, _ = distribution_sql(
                                                         'resistome',
                                                         'annotations',
                                                         'mutation_type',
                                                         mutation_type_threshold)

    t_keys = []
    t_values = []

    for key, val in zip(mutation_keys, mutation_values):
        if 'LARGE' not in key:
            t_keys.append(key)
            t_values.append(val)

    keys = t_keys
    mutation_values = t_values

    mutation_keys = convert_full_identifiers_to_abbreviations('mutations', keys)

    return journal_keys, journal_values, method_keys, method_values, mutation_keys, mutation_values


cur_obj = ConnectionManager(username=constants.DB_USERNAME,
                            password=constants.DB_PASSWORD,
                            db_name=constants.DB_NAME)
cursor = cur_obj.cursor