# enzymes.col: describes enzyme:reaction relationships
# reactions.dat can also be used in a similar way (probably better to use)
# genes dat: unique IDs for gene products + their synonyms
# proteins dat: unique ids for reactions from gene products
#
# parsing pattern:
#
# Key _ Value (note that value contains the delimiter '_' on many occasions)
#
# parsing rule:
#
# str.split("_",1), where 1 indicates at most 1 split at the first _ occurrence

# record separator: //

import json
import logging
import math
import os.path
import re
from collections import defaultdict

import psycopg2.extras

import resistome.sql.sql_build_utils as build_helper
from resistome import constants
from resistome.sql.sql_build_utils import SPECIES_SCHEMA, prepare_sql_query, insert_genome_sequences, \
    insert_sequence_data, parse_uniprot, insert_uniprot_data, insert_regulon_db_features, \
    parse_essential_genes_file, build_mutational_prediction_table, SNAP2_GENES, parse_fasta
from resistome.utils.regulondb_parser import parse_regulondb_distribution

log = logging.getLogger(__name__)

gene_data = 'genes.dat'
protein_data = 'proteins.dat'
rna_data = 'rnas.dat'
classes = 'classes.dat'
protein_sequence = 'protseq.fsa'
dna_sequence = 'dnaseq.fsa'
genome = 'genome.nt'
uniprot = 'uniprot-proteome_UP000000625.gff'
regulon_db = 'regulon_db_standardized.txt'
ecolinet = 'EcoliNet.v1.txt'
dna_binding_sites = 'dnabindsites.dat'
promoters = 'promoters.dat'
terminators = 'terminators.dat'
transcriptional_units = 'transunits.dat'
regulation = 'regulation.dat'


def biocyc_object_parser(obj, kv_delimiter, general_entry_delimiter, specific_entry_delimiter):

    """
    
    Parses the text flat record files that Biocyc uses to encode their data.
    
    :param obj: list, lines comprising the object
    :param kv_delimiter: str, denoting the separator between key/value
    :param general_entry_delimiter: str, default delimiter for separating keys with multiple values
    :param specific_entry_delimiter: str, override for general_entry_delimiter
    :return: dict (str : str) parsed object
    """

    obj_dict = {}
    previous_key = ''

    index = 0

    log.info('Parsing file into list of dicts (biocyc_object_parser)')

    for line in obj:

        # this leaves any empty entries intact
        line = line.translate(str.maketrans({'\n': None, '\r': None}))

        # removes html tags from input database
        line = re.sub('<.*?>', '', line)

        # continuation from previous key ('/' means line break in the file)
        if '/' == line[0]:

            try:
                obj_dict[previous_key] = obj_dict[previous_key] + \
                                                general_entry_delimiter + line.translate(str.maketrans({'/': None}))
            except:
                log.error('Error processing multi-line comment: %s (biocyc_object_parser)' % line)
                raise AssertionError('Error processing multi-line comment: %s' % line)

        else:
            # split only at first hyphen
            tokens = line.split(kv_delimiter, 1)

            try:
                # originally I was using biocyc field names as columns, postgres doesn't allow hyphens in column names.
                # you can remove this replace statement if desired.
                key = tokens[0].strip().replace("-", "_")
                value = tokens[1].strip()
            except:
                log.error('Error processing multi-line comment: %s (biocyc_object_parser)' % line)
                raise AssertionError('Failed to process single line comment: %s' % line)

            # multiple entries for same key
            if key in obj_dict:
                old_value = obj_dict[key]

                if key in specific_entry_delimiter:
                    new_value = old_value + specific_entry_delimiter[key] + value
                else:
                    new_value = old_value + general_entry_delimiter + value

                obj_dict[key] = new_value
            else:
                obj_dict[key] = value

            # reaction database has coefficients stored oddly
            # will result in a list of coefficients matching the react order
            if key == 'LEFT' or key == 'RIGHT':
                # still more reactants/products to go
                if index + 1 < len(obj) and obj[index + 1].find('^COEFFICIENT') == -1:
                    obj.insert(index + 1, '^COEFFICIENT - 1')
                # last entry in reaction database
                if index + 1 >= len(obj):
                    obj.append('^COEFFICIENT - 1')

            previous_key = key

        index = index + 1

    log.info('Finished parsing list into list of dicts (biocyc_object_parser)')

    return obj_dict


def process_biocyc_file(file_name):

    """
    
    Handles the parsing of biocyc flat files into list of lists that are subsequently parsed into dicts.
    
    :param file_name: 
    :return: 
    """

    fhandle = open(file_name, 'r')
    lines = fhandle.readlines()
    fhandle.close()

    log.debug('Loaded %s biocyc file for parsing (process_biocyc_file)' % file_name)

    biocyc_objects = []
    object_wrapper = []

    for line in lines:

        # ignore comments in biocyc flat files
        if line[0] != '#' and len(line) > 0:

            # record element, add to objectWrapper array so long as we haven't found the record separator
            if line[0:2] != "//":
                object_wrapper.append(line)

            # reached record delimiter
            else:
                biocyc_objects.append(object_wrapper)
                object_wrapper = []

    if len(object_wrapper) > 0:
        biocyc_objects.append(object_wrapper)

    return biocyc_objects


def extract_uniprot_xrefs(xref_text):

    """
    
    Extracts Uniprot cross references from dblinks field in biocyc entries and returns them as a list.
    
    :param xref_text: 
    :return: 
    """

    if 'UNIPROT' not in xref_text:

        return []

    # remove parantheses
    xref_text = xref_text.replace('(', '').replace(')', '')

    xrefs = xref_text.split(',')

    temp_output = []

    for token in xrefs:
        components = token.split(' ')
        database = components[0]
        accession = components[1].strip('"')
        temp_output.append((database, accession))

    # remove any non-uniprot cross-links
    temp_output = list(filter(lambda x: x[0] == 'UNIPROT', temp_output))

    return temp_output


def extract_regulator_data(regulatory_objects):

    """
    
    Extracts gene and other element transcriptional regulation from biocyc.
    
    :param regulatory_objects: 
    :return: 
    """

    output = []

    for regulatory_obj in regulatory_objects:

        mode = regulatory_obj.get('MODE', 'U')
        regulated_entity = regulatory_obj['REGULATED_ENTITY']
        if 'REGULATOR' in regulatory_obj:
            regulator = regulatory_obj['REGULATOR']
        else:
            print('skipped because it is missing the regulator entry', regulatory_obj)
            continue

        output.append((mode, regulated_entity, regulator))

    return output


def build_protein_disambiguator(protein_objects):

    """
    
    Protein objects are annotated in many different ways in Biocyc databases; some as the monomeric proteins,
    some as protein complexes, some as proteins with post-translational modifications, etc. This method attempts
    to convert an arbitrary protein id back to its monomeric equivalent(s).
    
    :param protein_objects: list of protein IDs
    :return: list of monomeric protein ids
    """

    name_converter = defaultdict(list)
    complexes = []

    monomeric_proteins = set()

    for unique_id in protein_objects:

        if 'UNMODIFIED_FORM' in protein_objects[unique_id] and 'Complexes' not in protein_objects[unique_id]['TYPES']:
            if 'CPLX' in protein_objects[unique_id]['UNMODIFIED_FORM']:
                complexes.append((unique_id, unique_id))
            else:
                name_converter[unique_id].append(protein_objects[unique_id]['UNMODIFIED_FORM'])
        elif 'Complexes' in protein_objects[unique_id]['TYPES'] or 'COMPONENTS' in protein_objects[unique_id]:
            complexes.append((unique_id, unique_id))
        else:
            name_converter[unique_id].append(unique_id)

    monomeric_proteins.update(name_converter.keys())

    seen = set()

    # this is like a low rent stack
    for (original_id, unique_id) in complexes:

        if unique_id not in protein_objects:
            continue

        if (original_id, unique_id) in seen:
            # this happens due to some weird complex annotations (self-referential complexes with no components?)
            continue

        seen.add((original_id, unique_id))
        components = set(protein_objects[unique_id].get('COMPONENTS', unique_id).split(' '))

        if 'UNMODIFIED_FORM' in protein_objects[unique_id]:
            components.add(protein_objects[unique_id]['UNMODIFIED_FORM'])

        monomeric_components = set(filter(lambda x: x in monomeric_proteins, components))

        complex_components = components - monomeric_components
        complex_components = [(original_id, x) for x in complex_components]

        complexes.extend(complex_components)

        for name in monomeric_components:
            name_converter[original_id].append(name)

    return name_converter


def parse_biocyc_flatentry(obj, file_type, kv_delimiter, entry_delimiter):

    """
    
    Wrapper for calling biocyc_object_parser. Automatically specifies some additional delimiters to better
    handle specific fields (see below for examples).
    
    :param obj: list of text representing a record
    :param file_type: str, name of the file the data is extracted from {'reactions', 'compounds'} have special handling
    :param kv_delimiter: str, key value delimiter
    :param entry_delimiter: str, delimiter for keys with multiple values
    :return: dict(str:str), parsed object
    """

    delim_map = {}

    if file_type in {'reactions'}:
        delim_map['^COEFFICIENT'] = ','
        delim_map['LEFT'] = '+'
        delim_map['RIGHT'] = '+'
        delim_map['IN_PATHWAY'] = ','
        delim_map['ENZYMATIC_REACTION'] = ','
    elif file_type in {'compounds'}:
        delim_map['CHEMICAL_FORMULA'] = ','
        delim_map['SYNONYMS'] = ','
        delim_map['REGULATES'] = ','
        delim_map['COFACTORS_OF'] = ','
    else:
        delim_map['SYNONYMS'] = ','
        delim_map['DBLINKS'] = ','
        delim_map['GO_TERMS '] = ','

    obj_dict = biocyc_object_parser(obj, kv_delimiter, entry_delimiter,
                                    delim_map)

    return obj_dict


def biocyc_name_to_gene_feature(object_dict, feature_name):

    """
    
    Attempts to convert the provided feature name into the associated genes, at the cost of checking the entire
    object dict christmas tree hierarchy.
    
    Searches: proteins, genes, promoters, terminators, and transcriptional units.
    
    :param object_dict: dict (str : dict (str : str)) representing all the data from a biocyc file
    :param feature_name: str, unique id of the feature in question
    :return: list of gene names (if found, otherwise empty list)
    """

    genes_found = []
    transcription_units_to_convert = []

    if feature_name in object_dict['proteins']:
        if 'GENE' in object_dict['proteins'][feature_name]:
            genes_found.append(object_dict['proteins'][feature_name]['GENE'])

    if feature_name in object_dict['genes']:
        genes_found.append(feature_name)

    if feature_name in object_dict['promoters']:

        if 'COMPONENT_OF' in object_dict['promoters'][feature_name]:
            transcription_units_to_convert.extend(object_dict['promoters'][feature_name]['COMPONENT_OF'].split(' '))

        # the remainder seem to be computationally inferred but unverified promoters that I will ignore.
        # especially since no location is given.

    if feature_name in object_dict['terminators']:

        if 'COMPONENT_OF' in object_dict['terminators'][feature_name]:
            transcription_units_to_convert.extend(object_dict['terminators'][feature_name]['COMPONENT_OF'].split(' '))

    if feature_name in object_dict['tu']:
        transcription_units_to_convert.append(feature_name)

    for tu in transcription_units_to_convert:

        if tu in object_dict['tu']:
            components = object_dict['tu'][tu]['COMPONENTS'].split(' ')

            # filter out non-gene elements
            genes = [x for x in components if x in object_dict['genes']]
            genes_found.extend(genes)

    # convert genes to accessions (preferred text id in database)

    converted_gene_names = []

    for gene in genes_found:

        if gene in object_dict['genes'] and 'ACCESSION_1' in object_dict['genes'][gene]:
            converted_gene_names.append(object_dict['genes'][gene].get('ACCESSION_1', gene))

    return converted_gene_names


def build_schema_dicts(schema_file_dict):

    """
    
    Handles the loading of the required biocyc files to build a species database.
    
    :param schema_file_dict: dict (str : str), file name to file path converter
    :return: 
    """

    entry_delimiter = ' '
    kv_delimiter = ' - '

    schema_dict = defaultdict(dict)

    for schema in schema_file_dict:

        files = schema_file_dict[schema]

        for (file_name, data_type) in files:

            if '.dat' in file_name:

                try:
                    biocyc_objects = process_biocyc_file(file_name)
                    biocyc_dict = dict()
                    for obj in biocyc_objects:
                        obj_dict = parse_biocyc_flatentry(obj, data_type, kv_delimiter, entry_delimiter)
                        biocyc_dict[obj_dict['UNIQUE_ID']] = obj_dict
                except IOError:
                    biocyc_dict = dict()

                schema_dict[schema][data_type] = biocyc_dict

            elif data_type in {'aa', 'dna', 'genome'}:

                schema_dict[schema][data_type] = parse_fasta(file_name)

            elif data_type in {'uniprot'}:

                try:
                    uniprot_dict = parse_uniprot(file_name)
                except IOError:
                    uniprot_dict = dict()

                schema_dict[schema][data_type] = uniprot_dict

            elif data_type in {'ecolinet_interactions'}:

                schema_dict[schema][data_type] = parse_ecolinet(file_name)

            elif data_type in {'regulondb_interactions'}:

                schema_dict[schema][data_type] = parse_regulondb_distribution()

    return schema_dict


def construct_schema_path_dict(source_dirs):

    """
    
    Sets up a dict describing where all the critical files are for processing a species database.
    
    :param source_dirs: 
    :return: 
    """

    schema_files = defaultdict(list)

    for (source, schema) in source_dirs:

        schema_files[schema].append((os.path.join(source, gene_data), 'genes'))
        schema_files[schema].append((os.path.join(source, protein_data), 'proteins'))
        schema_files[schema].append((os.path.join(source, rna_data), 'rna'))
        schema_files[schema].append((os.path.join(source, dna_sequence), 'dna'))
        schema_files[schema].append((os.path.join(source, protein_sequence), 'aa'))

        if os.path.exists(os.path.join(source, genome)):
            schema_files[schema].append((os.path.join(source, genome), 'genome'))
        elif os.path.exists(os.path.join(source, 'ecobase.nt')):
            schema_files[schema].append((os.path.join(source, 'ecobase.nt'), 'genome'))
        else:
            raise AssertionError('Did not find genome file for Biocyc distribution! %s' % source)

        schema_files[schema].append((os.path.join(source, classes), 'classes'))
        schema_files[schema].append((os.path.join(source, transcriptional_units), 'tu'))
        schema_files[schema].append((os.path.join(source, dna_binding_sites), 'dna_binding_sites'))
        schema_files[schema].append((os.path.join(source, promoters), 'promoters'))
        schema_files[schema].append((os.path.join(source, terminators), 'terminators'))
        schema_files[schema].append((os.path.join(source, regulation), 'regulation'))

        if schema == 'mg1655':
            schema_files[schema].append((os.path.join(constants.INPUT_DIR, 'biocyc', 'generic', uniprot), 'uniprot'))
            schema_files[schema].append((os.path.join(constants.INPUT_DIR, 'interactions', regulon_db),
                                         'regulondb_interactions'))
            schema_files[schema].append((os.path.join(constants.INPUT_DIR, 'interactions', ecolinet),
                                         'ecolinet_interactions'))

    return schema_files


def insert_genomic_features(cursor, schema, schema_entries, schema_objects,
                            strain_id, accession_to_gene_id, protein_name_converter):

    """
    
    Extracts the following genomic features:
    
    dnabindsites (dna binding sites, can get regulated genes from regulation.dat/transunits.dat
    promoters
    terminators
    
    Schema (per species)
    
    primary key
    unique id
    name
    feature type (terminator, promoter, binding site)
    start
    stop
    
    linked regulated by table (primary key, unique id, regulator)
    linked controlled genes (primary key, unique id, regulated gene)
    
    :param cursor: 
    :param schema: 
    :param schema_entries: 
    :param schema_objects: 
    :param protein_name_converter: 
    :return: 
    """

    objects_to_insert = []

    for unique_id in schema_objects['dna_binding_sites']:

        object_dict = schema_objects['dna_binding_sites'][unique_id]
        bs_type = 'dna_binding_site'
        if 'INVOLVED_IN_REGULATION' in object_dict:
            regulatory_references = object_dict['INVOLVED_IN_REGULATION'].split(' ')
        else:
            regulatory_references = []

        if 'LEFT_END_POSITION' in object_dict:
            start = int(object_dict['LEFT_END_POSITION'])
            stop = int(object_dict['RIGHT_END_POSITION'])
        elif 'ABS_CENTER_POS' in object_dict:
            try:
                center_position = float(object_dict['ABS_CENTER_POS'])
                site_length = float(object_dict.get('SITE_LENGTH', 1))
                start = int(math.floor(center_position - site_length / 2.0))
                stop = int(math.ceil(center_position - site_length / 2.0))
            except ValueError:
                # unknown position or position not given?
                start = -1
                stop = -1
        else:
            # if no known location is given
            start = -1
            stop = -1

        regulated_info = extract_regulator_data([schema_objects['regulation'][x] for x in regulatory_references])

        for (mode, regulated_entities, regulator) in regulated_info:

            converted_regulated_names = protein_name_converter.get(regulated_entities, [regulated_entities])
            converted_regulator_names = protein_name_converter.get(regulator, [regulator])

            regulators = []
            regulated = []

            for regulated_entity in converted_regulated_names:
                regulated.extend(biocyc_name_to_gene_feature(schema_objects, regulated_entity))

            for regulator_entity in converted_regulator_names:
                regulators.extend(biocyc_name_to_gene_feature(schema_objects, regulator_entity))

            objects_to_insert.append((unique_id, bs_type, start, stop, object_dict, regulated, regulators))

    for unique_id in schema_objects['promoters']:

        object_dict = schema_objects['promoters'][unique_id]
        site_type = 'promoter'
        if 'REGULATED_BY' in object_dict:
            regulatory_references = object_dict['REGULATED_BY'].split(' ')
        else:
            regulatory_references = []

        # transcriptional start site
        start = int(object_dict.get('ABSOLUTE_PLUS_1_POS', -1))
        stop = start + 1

        regulated_info = extract_regulator_data([schema_objects['regulation'][x] for x in regulatory_references])

        for (mode, regulated_entities, regulator) in regulated_info:

            converted_regulated_names = protein_name_converter.get(regulated_entities, [regulated_entities])
            converted_regulator_names = protein_name_converter.get(regulator, [regulator])

            regulators = []
            regulated = []

            for regulated_entity in converted_regulated_names:
                regulated.extend(biocyc_name_to_gene_feature(schema_objects, regulated_entity))

            for regulator_entity in converted_regulator_names:
                regulators.extend(biocyc_name_to_gene_feature(schema_objects, regulator_entity))

            objects_to_insert.append((unique_id, site_type, start, stop, object_dict, regulated, regulators))

    for unique_id in schema_objects['terminators']:

        object_dict = schema_objects['terminators'][unique_id]
        site_type = 'terminator'
        start = int(object_dict['LEFT_END_POSITION'])
        stop = int(object_dict['RIGHT_END_POSITION'])

        if 'Rho-Dependent-Terminators' in object_dict['TYPES']:
            # rho accession in mg1655
            regulator = ['B3783']
        else:
            regulator = []

        components = [x for x in object_dict['COMPONENT_OF'].split(' ') if 'TU' in x]

        gene_names = []

        for component in components:
            gene_names.extend(biocyc_name_to_gene_feature(schema_objects, component))

        objects_to_insert.append((unique_id, site_type, start, stop, object_dict, gene_names, regulator))

    already_inserted = dict()

    for (uid, site_type, start, stop, associated_data, regulated_genes, regulators) in objects_to_insert:

        if uid not in already_inserted:

            sql = prepare_sql_query('genomic_features',
                                    schema,
                                    schema_entries['genomic_features']) + ' RETURNING feature_id'

            cursor.execute(sql, (strain_id, uid, site_type, start, stop, json.dumps(associated_data), 'biocyc'))
            already_inserted[uid] = cursor.fetchone()['feature_id']

        feature_id = already_inserted[uid]

        sql = prepare_sql_query('genomic_feature_association', schema, schema_entries['genomic_feature_association'])

        for regulator in regulators:
            cursor.execute(sql, (accession_to_gene_id.get(regulator.upper(), None),
                                 feature_id, 'regulator', regulator.upper()))

        for regulated_gene in regulated_genes:
            cursor.execute(sql, (accession_to_gene_id.get(regulated_gene.upper(), None),
                                 feature_id, 'regulated', regulated_gene.upper()))


def extract_go_information(schema, strain_objects, root_nodes, seen_go_codes):

    go_terms_to_include = set()

    for class_key in strain_objects[schema]['classes']:

        if 'GO:' not in class_key or class_key in seen_go_codes:
            continue

        seen_go_codes.add(class_key)

        types = strain_objects[schema]['classes'][class_key]['TYPES'].split(' ')

        if class_key in root_nodes:
            types = []

        name = strain_objects[schema]['classes'][class_key]['COMMON_NAME']

        if 'Gene-Ontology-Terms' not in types:
            go_terms_to_include.add(class_key)

        for ancestor in types:
            if 'GO:' not in ancestor:
                continue
            go_terms_to_include.add(ancestor)

    return go_terms_to_include, seen_go_codes


def parse_ecolinet(filename):
    """

    Parses the co-functional relation interaction network developed by Kim et al. (10.1093/database/bav001).

    :param filename:
    :return:
    """

    output = set()

    with open(filename, 'r') as f:
        for line in f:
            tokens = line.strip().split('\t')
            geneA = tokens[0].upper()
            geneB = tokens[1].upper()
            direction = '?'

            output.add(('protein-protein', geneA, geneB, direction, 'ecolinet_v1'))

    return output


def insert_interaction_data(cur, schema, table, columns, strain_id, data_tuples):

    for data_tuple in data_tuples:
        sql = prepare_sql_query(table,
                                schema,
                                columns,
                                ['%s'] * len(SPECIES_SCHEMA['interactions']))
        cur.execute(sql, (strain_id, ) + data_tuple)


def main(database_cursor: psycopg2._psycopg.cursor,
         biocyc_directory: str,
         species: str,
         strain_of_interest: str):

    """
    
    Driver method for database construction; does the heavy lifting of getting the data sources together and 
    then drives the inserts.
    
    :return: 
    """

    raise NotImplementedError('Deprecated. You should try to rely as much as possible on public annotation sources '
                              'that can be used freely (e.g. NCBI).')


    strain_id = build_helper.insert_strain(database_cursor, species, strain_of_interest)

    manual_gene_locations = defaultdict(list)

    # required for using Biocyc 24.1 MG1655 database since many of the pseudo genes are referenced elsewhere
    with open(os.path.join(constants.INPUT_DIR, 'biological_info', 'pseudogene_locations.txt')) as f:
        f.readline()
        for line in f:
            tokens = line.strip().split('\t')
            location = tokens[-1].split('..')
            manual_gene_locations[tokens[0]].append(tuple([int(x) for x in location]))

    schema_files = construct_schema_path_dict([(biocyc_directory, strain_of_interest)])

    converter = {}
    # name in column names (below) converted to the biocyc equivalent
    converter['COEFFICIENT'] = '^COEFFICIENT'
    converter['RIGHT_HS'] = 'RIGHT'
    converter['LEFT_HS'] = 'LEFT'
    converter['ORPHAN'] = 'ORPHAN?'

    strain_objects = build_schema_dicts(schema_files)

    # go data
    seen_go_codes = set()
    root_nodes = {'GO:0008150', 'GO:0005575', 'GO:0003674'}
    go_terms_to_include = set()

    # generates schema
    unique_id_to_accession = dict()

    # load essential genes
    essential_set = parse_essential_genes_file(os.path.join(constants.INPUT_DIR,
                                                            'biological_info',
                                                            'essential_genes_ecoli.txt'))

    accession_to_gene_id = dict()

    # get existing GO ids
    database_cursor.execute('select go_id, go_term from go_table')

    go_term_to_go_id = dict()
    for record in database_cursor:
        go_term_to_go_id[record['go_term']] = record['go_id']

    for schema in strain_objects:

        go_terms_extracted, seen_updated = extract_go_information(schema,
                                                                  strain_objects,
                                                                  root_nodes,
                                                                  seen_go_codes)
        go_terms_to_include.update(go_terms_extracted)
        seen_go_codes.update(seen_updated)

        uniprot_to_accession = dict()
        gene_objects = strain_objects[schema]['genes'].values()

        protein_name_converter = build_protein_disambiguator(strain_objects[schema]['proteins'])

        for obj_dict in gene_objects:

            accession = obj_dict.get('ACCESSION_1', None)

            if accession is None:
                continue

            accession = accession.upper()

            unique_gene_id = obj_dict['UNIQUE_ID']
            start = obj_dict.get('LEFT_END_POSITION', None)
            stop = obj_dict.get('RIGHT_END_POSITION', None)
            direction = obj_dict.get('TRANSCRIPTION_DIRECTION', None)
            unique_id_to_accession[unique_gene_id] = (accession, start, stop, direction)

            common_name = obj_dict.get('COMMON_NAME', unique_gene_id)
            essential = common_name.upper() in essential_set

            synonyms = []

            if 'SYNONYMS' in obj_dict:
                synonyms = obj_dict['SYNONYMS'].split(',')

            if 'ACCESSION_2' in obj_dict:
                synonyms.append(obj_dict['ACCESSION_2'])

            protein_ids = []

            uniprot_refs = []

            if 'PRODUCT' in obj_dict:
                protein_ids.extend(obj_dict['PRODUCT'].replace('|', '').split(' '))
                for id in protein_ids:
                    unique_id_to_accession[id] = (accession, None, None, None)
                    if id in strain_objects[schema]['proteins']:
                        uniprot_refs.extend(extract_uniprot_xrefs(strain_objects[schema]['proteins'][id].get('DBLINKS', '')))
                    elif id in strain_objects[schema]['rna']:
                        uniprot_refs.extend(extract_uniprot_xrefs(strain_objects[schema]['rna'][id].get('DBLINKS', '')))

            if len(uniprot_refs) > 0:
                for (_, uniprot_id) in uniprot_refs:
                    uniprot_to_accession[uniprot_id] = accession

            go_terms = []

            for protein_id in protein_ids:

                if protein_id in strain_objects[schema]['rna']:
                    table = 'rna'
                    protein_object_dict = strain_objects[schema][table][protein_id]
                    if 'GO_TERMS' in protein_object_dict:
                        go_terms.extend(protein_object_dict['GO_TERMS'].replace('|', '').split(' '))
                elif protein_id in strain_objects[schema]['proteins']:
                    table = 'proteins'
                    protein_object_dict = strain_objects[schema][table][protein_id]
                    if 'GO_TERMS' in protein_object_dict:
                        go_terms.extend(protein_object_dict['GO_TERMS'].replace('|', '').split(' '))

            # start - 1 converts the index to zero

            arguments = [strain_id,
                         unique_gene_id,
                         common_name,
                         synonyms,
                         protein_ids,
                         accession,
                         # pseudo gene status
                         True if start is None and stop is None else False,
                         essential]

            SQL = prepare_sql_query('genes', schema, SPECIES_SCHEMA['genes'], arguments)
            SQL += ' RETURNING gene_id'

            if start is not None and stop is not None:
                database_cursor.execute(SQL, arguments)
                gene_id = database_cursor.fetchone()['gene_id']
                database_cursor.execute('INSERT INTO gene_locations (gene_id, start, stop, direction) '
                                        'VALUES (%s, %s, %s, %s)',
                            (gene_id, start, stop, direction))
            elif unique_gene_id in manual_gene_locations:
                print('Gene %s appears to be a pseudogene; using manual gene locations to insert' % unique_gene_id)
                locations = manual_gene_locations[unique_gene_id]
                database_cursor.execute(SQL, arguments)
                gene_id = database_cursor.fetchone()['gene_id']
                for (start_x, stop_y) in locations:
                    database_cursor.execute('INSERT INTO gene_locations (gene_id, start, stop, direction) VALUES (%s, %s, %s, %s)',
                                            (gene_id, start_x, stop_y, direction))
            else:
                print("Missing gene boundaries? Might be a pseudogene?", arguments)
                continue

            accession_to_gene_id[accession] = gene_id

            for go_term in go_terms:

                if go_term not in go_terms_to_include:
                    continue

                SQL = prepare_sql_query('go_terms', schema, SPECIES_SCHEMA['go_terms'], [gene_id, go_term])
                database_cursor.execute(SQL, (gene_id, go_term_to_go_id[go_term]))

            if unique_gene_id in strain_objects[schema]['dna']:
                insert_sequence_data(database_cursor,
                                     schema,
                                     SPECIES_SCHEMA['dna'],
                                     'nt_sequences',
                                     gene_id,
                                     strain_objects[schema]['dna'][unique_gene_id])

            for protein_id in protein_ids:
                if protein_id in strain_objects[schema]['aa']:
                    insert_sequence_data(database_cursor,
                                         schema,
                                         SPECIES_SCHEMA['aa'],
                                         'aa_sequences',
                                         gene_id,
                                         strain_objects[schema]['aa'][protein_id])

        insert_genome_sequences(database_cursor, schema,
                                SPECIES_SCHEMA['genome'],
                                strain_id,
                                strain_objects[schema]['genome'])

        insert_interaction_data(database_cursor, schema,
                                'interactions',
                                SPECIES_SCHEMA['interactions'],
                                strain_id,
                                strain_objects[schema].get('ecolinet_interactions', []))

        insert_genomic_features(database_cursor,
                                schema, SPECIES_SCHEMA, strain_objects[schema],
                                strain_id,
                                accession_to_gene_id,
                                protein_name_converter)

        if 'uniprot' in strain_objects[schema]:
            uniprot_dict = strain_objects[schema]['uniprot']
            insert_uniprot_data(database_cursor, schema, SPECIES_SCHEMA['uniprot'], uniprot_dict,
                                uniprot_to_accession,
                                accession_to_gene_id)

        insert_regulon_db_features(accession_to_gene_id, database_cursor, schema, strain_id, strain_objects)

    # allowable method types include 'inps' and 'snap2'
    # warning: loading the snap2 data takes quite a bit of time...
    # note: I am manually encoding which genes to add to the db. you can change the method to just load everything
    # but the problem is that this takes a huge amount of space and usually is not necessary.
    build_mutational_prediction_table(cur=database_cursor,
                                      strain_to_process=strain_of_interest,
                                      unique_id_to_accession_dict=unique_id_to_accession,
                                      accession_to_gene_id=accession_to_gene_id,
                                      methods={'inps', 'snap2'},
                                      variant_effect_predictor_genes=SNAP2_GENES)

    print('Finished parsing/uploading provided biocyc database: %s' % strain_of_interest)


if __name__ == '__main__':

    # todo: add table of dna binding proteins - consensus sequences

    # manually get cursor to avoid depending on Database utils for this
    # also, you need to run the biocyc_data_parser to build the biocyc tables first

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    root_dir = constants.INPUT_DIR
    ECOLI_DIR = os.path.join(root_dir, 'biocyc', 'mg1655', 'data')
    ECOLIB_DIR = os.path.join(root_dir, 'biocyc', 'rel606', 'data')
    ECOLIW_DIR = os.path.join(root_dir, 'biocyc', 'ecoliw', 'data')

    source_dirs_x = [(ECOLI_DIR, 'mg1655'),
                   (ECOLIB_DIR, 'rel606'),
                   (ECOLIW_DIR, 'w')]

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:

        with open(os.path.join(constants.INPUT_DIR, 'sql', 'resistome_biocyc_schema.sql'), 'r') as f:
            sql_schema = ''.join(f.readlines())
            # make schema species specific
            cursor.execute(sql_schema)

        for (data_dir, strain) in source_dirs_x:
            main(cursor, biocyc_directory=data_dir, species='Escherichia coli', strain_of_interest=strain)

    connect.commit()
    connect.close()



