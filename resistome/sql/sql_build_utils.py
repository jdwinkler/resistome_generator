import json
import os.path
from collections import defaultdict

import psycopg2
import psycopg2.sql as psql
import psycopg2.extras

from resistome import constants
from resistome.utils import database_utils
from typing import List, Tuple
import glob

from resistome.utils.regulondb_parser import extract_mg1655_genome_features

SPECIES_SCHEMA = dict()

# schema for storing species data
SPECIES_SCHEMA['genes'] = ['strain_id',
                           'biocyc_id',
                           'name',
                           'synonyms',
                           'products',
                           'accession',
                           'description',
                           'pseudogene',
                           'essential']
SPECIES_SCHEMA['dna'] = ['gene_id', 'nt_seq']
SPECIES_SCHEMA['aa'] = ['gene_id', 'aa_seq']
SPECIES_SCHEMA['genome'] = ['strain_id', 'chromosome', 'sequence']
SPECIES_SCHEMA['uniprot'] = ['gene_id', 'region', 'start', 'stop', 'note']
SPECIES_SCHEMA['go_terms'] = ['gene_id', 'go_id']
SPECIES_SCHEMA['interactions'] = ['strain_id',
                                  'interaction_type',
                                  'regulator_name',
                                  'target_name',
                                  'regulator',
                                  'target',
                                  'direction',
                                  'db_source']
SPECIES_SCHEMA['metabolomics'] = ['gene_id', 'metabolite_id']
SPECIES_SCHEMA['genomic_features'] = ['strain_id', 'biocyc_id', 'feature_type', 'start', 'stop', 'associated_data', 'source']
SPECIES_SCHEMA['genomic_feature_association'] = ['gene_id', 'feature_id', 'relationship', 'name']

# the Resistome schema
RESISTOME_SCHEMA = {'papers': ['title', 'doi', 'year', 'research_group', 'journal',
                         'methods', 'score', 'reason', 'genome_reference', 'designs', 'comments'],
              'paper_tags': ['paper_id', 'tag'],
              'mutants': ['paper_id', 'name',
                          'species', 'original_strain', 'strain', 'oxygenation',
                          'medium', 'carbon_source', 'media_supplements', 'ph',
                          'vessel_type', 'culture_volume', 'vessel_volume',
                          'temperature', 'rotation', 'fold_improvement',
                          'initial_fitness', 'final_fitness', 'fitness_unit',
                          'comments'],
              'mutant_methods': ['mutant_id', 'method'],
                    'mutations': ['paper_id', 'mutant_id', 'name', 'species', 'strain', 'effects', 'original'],
                    'expressions': ['paper_id', 'mutant_id', 'study_id', 'name', 'species', 'strain', 'status', 'pvalue',
                              'fold_change'],
                    'expression_studies': ['paper_id', 'mutant_id', 'accession', 'measurement_method', 'statistical_method',
                                     'duration', 'growth_phase', 'stressor_level', 'stressor_units'],
                    'annotations': ['gene_id', 'mutation_type', 'annotation'],
                    'phenotypes': ['mutant_id', 'phenotype', 'phenotype_class', 'phenotype_type', 'ontology_root',
                             'resistance_level', 'resistance_unit'],
                    'gene_standardization': ['strain', 'species_accession', 'mg1655_accession'],
                    'gene_metadata': ['name', 'species', 'strain', 'species_accession', 'mg1655_accession'],
                    'term_explanation': ['term_type', 'internal_name', 'explanation'],
                    'phenotype_standardization': ['standard_name', 'phenotype_type', 'root_class',
                                            'specific_classes'],
                    'abbreviations': ['entry', 'entry_type', 'converted_entry'],
                    'gene_ontology': ['accession', 'go_term']}

SNAP2_GENES = {'B3357',
               'B3251',
               'B3988',
               'B3261',
               'B0462',
               'B3987',
               'B3067',
               'B4018',
               'B3783',
               'B4063',
               'B1503',
               'B0464',
               'B3938',
               'B0002',
               'B1530',
               'B3340',
               'B3650',
               'B0635',
               'B2231'}


def insert_strain(cursor: psycopg2._psycopg.cursor, species: str, strain: str):

    cursor.execute('INSERT INTO strain (species, strain) VALUES (%s, %s) RETURNING strain_id',
                   (species, strain))

    return cursor.fetchone()['strain_id']


def prepare_sql_query(table, schema, columns, arguments=None):

    """

    Generates an sql insert query based on the provided inputs.

    :param schema: name of the schema (str), usually resistome
    :param table: str, name of the table
    :param columns: iterable of strs, get from column dict
    :param arguments
    :return: formatted SQL query
    """

    if arguments is None:
        arguments = columns

    fields = '(' + ', '.join(list(columns)) + ')'

    if schema == 'resistome':
        # hack to deal with bad decisions.
        SQL = 'INSERT INTO ' + schema + '.' + table + ' ' + fields + \
              ' VALUES (' + ', '.join(['%s'] * len(arguments)) + ')'
    else:
        SQL = 'INSERT INTO ' + table + ' ' + fields + ' VALUES (' + ', '.join(['%s'] * len(arguments)) + ')'

    return SQL


def prepare_tuple_style_sql_query(table, schema, columns):

    fields = '(' + ', '.join(list(columns)) + ')'

    if schema == 'resistome':
        SQL = 'INSERT INTO resistome.' + table + ' ' + fields + ' VALUES %s'
    else:
        SQL = 'INSERT INTO ' + table + ' ' + fields + ' VALUES %s'

    return SQL


def insert_genome_sequences(cur, schema, columns, strain_id, genome_object):

    for key in genome_object:
        gen_seq = genome_object[key]
        sql = prepare_sql_query('genome', schema, columns, ['%s', '%s', '%s'])
        cur.execute(sql, (strain_id, key, gen_seq))


def insert_sequence_data(cur, schema, columns, sequence_table, gene_id, sequence):

    assert sequence_table == 'aa_sequences' or sequence_table == 'nt_sequences'

    sql = prepare_sql_query(sequence_table, schema, columns, ['%s', '%s'])
    cur.execute(sql, (gene_id, sequence))


def parse_uniprot(file_path):

    """

    Extract annotations from an uniprot style file downloaded from their website.

    :param file_path:
    :return:
    """

    output = defaultdict(list)
    with open(file_path, 'r') as fhandle:
        for line in fhandle:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            uniprot_id = tokens[0]
            entry_type = tokens[2]
            start = tokens[3]
            stop = tokens[4]

            info = tokens[-1].split(';')
            info = list(filter(lambda x: 'Note=' in x, info))

            if len(info) == 0:
                note = None
            else:
                note = info[0].split('=')[1].strip()
            output[uniprot_id].append((entry_type, start, stop, note))
        return output


def insert_uniprot_data(cur, schema, columns, uniprot_dict, uniprot_to_accession, accession_to_gene_id):

    for uniprot_id in uniprot_to_accession:

        if uniprot_id not in uniprot_dict:
            continue

        accession = uniprot_to_accession[uniprot_id]

        uniprot_tuples = []
        sql = prepare_tuple_style_sql_query('uniprot', schema, columns)
        for (region, start, stop, note) in uniprot_dict[uniprot_id]:

            try:
                uniprot_tuples.append((accession_to_gene_id[accession],
                                       region, start - 1 if isinstance(start, int) else start, stop, note))
            except:
                raise

        psycopg2.extras.execute_values(cur, sql, uniprot_tuples, page_size=2000)


def insert_regulon_db_features(accession_to_gene_id, database_cursor, schema, strain_id, strain_objects=None):

    regulon_db_genome_features, associations = extract_mg1655_genome_features(strain_id)

    SQL = prepare_tuple_style_sql_query('genomic_features', schema, SPECIES_SCHEMA['genomic_features'])

    psycopg2.extras.execute_values(database_cursor, SQL, regulon_db_genome_features, page_size=2000)

    for feature_class, association_entry in associations.items():

        features = [(x[1], x[5]) for x in filter(lambda x: x[2] == feature_class, regulon_db_genome_features)]
        database_cursor.execute('select biocyc_id, feature_id from genomic_features where feature_type = %s',
                                (feature_class,))

        unique_id_to_feature_id = dict()
        for record in database_cursor:
            biocyc_id = record['biocyc_id']
            feature_id = record['feature_id']
            unique_id_to_feature_id[biocyc_id] = feature_id

        feature_tuples = []
        for unique_id, data_json in features:
            genes = json.loads(data_json)[association_entry]
            if len(genes) == 0:
                continue
            for gene in genes:
                feature_tuples.append((accession_to_gene_id.get(gene.upper(), None),
                                       unique_id_to_feature_id[unique_id],
                                       associations[feature_class].lower(),
                                       gene))
        SQL = prepare_tuple_style_sql_query('genomic_feature_association', schema,
                                            SPECIES_SCHEMA['genomic_feature_association'])

        psycopg2.extras.execute_values(database_cursor, SQL, feature_tuples, page_size=2000)


def parse_essential_genes_file(filename):

    """

    Parses list of essential genes.

    :param filename:
    :return:
    """

    essential_set = set()

    with open(filename, 'r') as f:

        f.readline()
        for line in f:
            tokens = line.strip().split()
            gene_name = tokens[0].upper()
            essential_set.add(gene_name)

    return essential_set


def build_mutational_prediction_table(cur,
                                      strain_to_process,
                                      unique_id_to_accession_dict,
                                      accession_to_gene_id,
                                      methods={'inps'}, variant_effect_predictor_genes=set()):

    """

    Loads data predicting the effect of mutations on proteins. Currently only loads INPS data by default.

    INPS data obtained from here: 10.1093/bioinformatics/btv291

    :param cur
    :param unique_id_to_accession_dict
    :param methods
    :return:
    """

    columns = ['gene_id', 'position', 'wt_aa', 'mutant_aa', 'score', 'method']

    if not os.path.exists(os.path.join(constants.INPUT_DIR, 'protein_predictions', 'biocyc_to_accession_map.txt')):
        raise AssertionError('Missing mapping file-run inputs/protein_predictions/match_seq.py to generate.')

    biocyc_to_accession = dict()
    with open(os.path.join(constants.INPUT_DIR, 'protein_predictions', 'biocyc_to_accession_map.txt')) as f:
        for line in f:
            tokens = line.strip().split('\t')
            biocyc_id = tokens[1]
            # globally unique
            accession = tokens[2]
            biocyc_to_accession[biocyc_id] = accession

    if 'inps' in methods:

        inps_protein_data = os.path.join(constants.INPUT_DIR, 'protein_predictions', 'inps.pred.txt')

        if not os.path.exists(inps_protein_data):
            print('Protein stability data for INPS was not found, skipping requested upload...')
            return

        with open(inps_protein_data, 'r') as f:

            for line in f:
                tokens = line.strip().split('\t')

                gene = tokens[0]
                aa_change = tokens[2]

                wt_aa = aa_change[0]
                mut_aa = aa_change[-1]
                position = int(aa_change[1:-1]) - 1
                score = tokens[3]

                sql = prepare_sql_query('protein_stability', strain_to_process, columns)
                if gene in accession_to_gene_id:
                    cur.execute(sql, (accession_to_gene_id[gene], position, wt_aa, mut_aa, score, 'INPS'))

    if 'snap2' in methods:

        if len(variant_effect_predictor_genes) == 0:
            raise AssertionError('Given the size of the SNAP2 dataset, you need to specify which genes you'
                                 'want to extract data from.')

        demask_dir = os.path.join(constants.INPUT_DIR, 'protein_predictions', 'snap2')

        target_directory = os.path.join(demask_dir, strain_to_process.upper())

        if not os.path.exists(target_directory):
            print('Did not find SNAP2 directory for strain, skipping: %s' % strain_to_process)
            return
        else:

            print('Uploading SNAP2 data for %s' % strain_to_process)

            files = os.listdir(target_directory)
            files = list(filter(lambda x: 'snap2' not in x and 'fasta' in x, files))
            sql = prepare_tuple_style_sql_query('protein_stability', strain_to_process, columns)

            for fname in files:

                sequence_id = fname
                snap2_data = fname + '.snap2'

                fasta = parse_fasta(os.path.join(target_directory, sequence_id))

                try:

                    accession_in_file = list(fasta.keys())[0]

                    try:
                        accession = biocyc_to_accession[accession_in_file]
                    except KeyError:
                        print('Expected all biocyc IDs to be in the biocyc => accession mapping! %s is not.'
                              % accession_in_file)
                        raise

                    if accession not in variant_effect_predictor_genes:
                        continue

                    if accession not in accession_to_gene_id:
                        print('unknown gene accession', accession)
                        continue

                    demask_tuples = parse_snap2_file(accession_to_gene_id[accession],
                                                    accession,
                                                    os.path.join(target_directory, snap2_data),
                                                    score_threshold=0.0)

                    assert len(demask_tuples) > 0

                    psycopg2.extras.execute_values(cur,
                                                   sql,
                                                   demask_tuples,
                                                   page_size=10000)
                except Exception as e:
                    print('Error when processing %s: %s; skipping' % (snap2_data, e))
                    continue

    if 'DEMASK' in methods:

        if len(variant_effect_predictor_genes) == 0:
            raise AssertionError('Given the size of the DeMaSk dataset, you need to specify which genes you'
                                 'want to extract data from.')

        demask_dir = os.path.join(constants.INPUT_DIR, 'protein_predictions', 'demask')

        target_directory = os.path.join(demask_dir, strain_to_process.upper())

        if not os.path.exists(target_directory):
            print('Did not find DeMaSk directory for strain, skipping: %s' % strain_to_process)
            return
        else:

            print('Uploading DeMaSk data for %s' % strain_to_process)

            sql = prepare_tuple_style_sql_query('protein_stability', strain_to_process, columns)

            for fname in glob.glob(os.path.join(target_directory, '*.txt')):

                accession_in_file = os.path.splitext(os.path.split(fname)[1])[0]

                try:

                    try:
                        accession = unique_id_to_accession_dict[accession_in_file]
                    except KeyError:
                        print('Expected all biocyc IDs to be in the biocyc => accession mapping! %s is not.'
                              % accession_in_file)
                        raise

                    if accession not in variant_effect_predictor_genes:
                        continue

                    if accession not in accession_to_gene_id:
                        print('unknown gene accession', accession)
                        continue

                    demask_tuples = parse_demask_file(accession_to_gene_id[accession],
                                                     accession)

                    assert len(demask_tuples) > 0

                    psycopg2.extras.execute_values(cur,
                                                   sql,
                                                   demask_tuples,
                                                   page_size=10000)
                except Exception as e:
                    print('Error when processing %s: %s; skipping' % (accession_in_file, e))
                    continue


def parse_fasta(filename):

    """

    Simple fasta parser.

    :param filename:
    :return:
    """

    fasta_dict = defaultdict(list)

    fasta_naive_dict = database_utils.load_fasta(filename)

    for key in fasta_naive_dict:

        if '|' in key:
            unique_id = key.strip().split('|')[2].split(' ')[0]
        else:
            unique_id = key.strip().split(' ')[0][1:]

        fasta_dict[unique_id] = fasta_naive_dict[key]

    return fasta_dict


def parse_snap2_file(gene_id, accession, filename, score_threshold=0.0):

    """

    Parses SNAP2 protein mutation effect prediction files. These data are essentially matrices with the
    sequence as rows and a fixed order of AAs as columns, and then a score is assigned using a neural network
    based on mutational training data. Not used currently due to the size of the datasets, but may be useful
    in the future.

    See here for details about SNAP2: https://rostlab.org/owiki/index.php/Snap2

    :param accession
    :param filename:
    :param score_threshold
    :return:
    """

    output_tuples = []

    with open(filename, 'r') as f:

        for line in f:

            if '=>' not in line:
                continue

            tokens = line.strip().split('=>')

            aa_pos = tokens[0].strip()

            score = tokens[1].split('sum =')[1].strip()

            wt_aa = aa_pos[0]
            position = str(int(aa_pos[1:-1])-1)
            candidate_replacement_aa = aa_pos[-1]

            factor = -1 if 'Neutral' in tokens[1] else 1

            score = float(score) / 100.0 * factor

            if abs(score) < score_threshold:
                continue
            else:
                output_tuples.append((gene_id, position, wt_aa, candidate_replacement_aa, score, 'SNAP2'))

    return output_tuples


def parse_demask_file(gene_id, filename) -> List[Tuple[str, ...]]:

    output_tuples = []

    with open(filename) as f:
        header = {x: k for k, x in enumerate(f.readline().strip().split('\t'))}
        for line in f:
            tokens = line.strip().split('\t')

            position = tokens[header['position']]
            wt_aa = tokens[header['WT']]
            mut_aa = tokens[header['var']]
            score = tokens[header['score']]

            output_tuples.append((gene_id, position, wt_aa, mut_aa, score, 'DeMaSk'))

    return output_tuples
