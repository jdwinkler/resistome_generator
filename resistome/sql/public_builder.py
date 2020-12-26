import itertools
import os.path
from collections import defaultdict

import gzip
import networkx
import obonet
import psycopg2
import psycopg2.extras

from resistome import constants
from resistome.sql.sql_build_utils import prepare_sql_query, prepare_tuple_style_sql_query, SPECIES_SCHEMA
from resistome.utils import database_utils


def insert_genetic_code(cur, columns, filename):
    """

    Loads E. coli genetic code into table.

    :param cur:
    :param columns:
    :param filename:
    :return:
    """

    with open(filename, 'r') as fhandle:
        fhandle.readline()
        for line in fhandle:
            tokens = line.strip().split('\t')
            codon = tokens[0]
            aa = tokens[1]
            frequency = tokens[2]

            sql = prepare_sql_query('genetic_code', 'public', columns,
                                    ['%s'] * len(columns))
            cur.execute(sql, (11, codon, aa, frequency))


def insert_interaction_data(cur, schema, table, columns, data_tuples):

    for data_tuple in data_tuples:
        sql = prepare_sql_query(table,
                                schema,
                                columns,
                                ['%s'] * len(columns))
        cur.execute(sql, data_tuple)


def main(cur):

    # adds genetic code to public table
    insert_genetic_code(cur,
                        ['code', 'codon', 'aa', 'frequency'],
                        os.path.join(constants.INPUT_DIR,
                                     'biological_info',
                                     'Genetic Code EC.txt'))

    graph = obonet.read_obo(os.path.join(constants.INPUT_DIR, 'ontologies', 'go_DEC2020.obo'),
                            ignore_obsolete=False)

    adj_list = set()
    go_term_table = []

    alt_ids = dict()

    for node, data in graph.nodes(data=True):
        go_term_table.append((node, [x for x in networkx.ancestors(graph, node)], data['name']))
        if 'alt_id' in data:
            alt_id_names = data['alt_id']
            for alt_id in alt_id_names:
                go_term_table.append((alt_id, [x for x in networkx.ancestors(graph, node)], data['name']))
            alt_ids[node] = alt_id_names

    for (p, c, _) in graph.edges:
        adj_list.add((p, c))

        if p in alt_ids:
            for x in alt_ids[p]:
                adj_list.add((x, c))
        if c in alt_ids:
            for y in alt_ids[c]:
                adj_list.add((p, y))
        if p in alt_ids and c in alt_ids:
            for z, a in itertools.product(alt_ids[p], alt_ids[c]):
                adj_list.add((z, a))

    SQL = prepare_tuple_style_sql_query('go_table', 'public', ['go_term', 'ancestors', 'name']) + \
          ' returning go_id, go_term'
    results = psycopg2.extras.execute_values(cur, sql=SQL, argslist=go_term_table, page_size=2000,
                                             fetch=True)
    go_to_id = dict()
    for record in results:
        go_to_id[record['go_term']] = record['go_id']

    go_adj_with_pk = []
    for (p, c) in adj_list:
        go_adj_with_pk.append((p, c, go_to_id[p], go_to_id[c]))

    SQL = prepare_tuple_style_sql_query('go_adj_list', 'public', ['parent', 'child', 'parent_id', 'child_id'
                                                                                                  ''])
    psycopg2.extras.execute_values(cur, sql=SQL, argslist=go_adj_with_pk, page_size=2000)


def build_sauer_metabolite_table(cur, directory_path):

    """

    Loads Sauer gene deletion-metabolite data files from supplementary data. Reference: 10.15252/msb.20167150

    :param cur:
    :return:
    """

    # z-score minimum for inclusion in the table
    z_score_threshold = 2.765

    # Note: you wan to call this function after you finish all database updates.
    unified_std = database_utils.get_unified_standard()

    gene_annotations = []
    # order of data
    with open(os.path.join(directory_path, 'sample_id_modzscore.txt'), 'r') as f:
        for line in f:
            gene_annotations.append(line.strip())

    neg_ion_annotations = []
    pos_ion_annotations = []

    sql = prepare_tuple_style_sql_query('ms_ions', 'public', ['ion_id', 'kegg_id', 'name']) + ' returning ion_pk, ion_id'

    # metabolite information
    with open(os.path.join(directory_path, 'annotation_pos.txt'), 'r') as f:

        for line in f:
            tokens = line.strip().split('\t')
            ion_mass = tokens[0]
            if len(tokens) > 2:
                name = tokens[1]
                kegg_id = tokens[2]
            else:
                name = None
                kegg_id = None
            pos_ion_annotations.append((ion_mass, kegg_id, name))
            # cur.execute(sql, (ion_mass, kegg_id, name))

    with open(os.path.join(directory_path, 'annotation_neg.txt'), 'r') as f:

        for line in f:
            tokens = line.strip().split('\t')
            ion_mass = tokens[0]
            if len(tokens) > 2:
                name = tokens[1]
                kegg_id = tokens[2]
            else:
                name = None
                kegg_id = None
            neg_ion_annotations.append((ion_mass, kegg_id, name))
            # cur.execute(sql, (ion_mass, kegg_id, name))

    ion_mass_to_pk = dict()
    results = psycopg2.extras.execute_values(cur, sql, argslist=neg_ion_annotations + pos_ion_annotations,
                                             page_size=2000,
                                             fetch=True)

    for result in results:
        ion_mass_to_pk[result['ion_id']] = result['ion_pk']

    neg_matrix = []
    pos_matrix = []

    with gzip.open(os.path.join(directory_path, 'modzscore_neg.tsv.gz'), 'rt') as f:
        for line in f:
            metabolite_zscores = list(map(float, line.strip().split('\t')))
            neg_matrix.append(metabolite_zscores)

    with gzip.open(os.path.join(directory_path, 'modzscore_pos.tsv.gz'), 'rt') as f:
        for line in f:
            metabolite_zscores = list(map(float, line.strip().split('\t')))
            pos_matrix.append(metabolite_zscores)

    pos_affected_metabolites = defaultdict(list)
    neg_affected_metabolites = defaultdict(list)

    for j in range(0, len(gene_annotations)):
        for i in range(0, len(neg_ion_annotations)):
            if abs(neg_matrix[i][j]) > z_score_threshold:
                neg_affected_metabolites[gene_annotations[j]].append((i, neg_matrix[i][j]))

    for j in range(0, len(gene_annotations)):
        for i in range(0, len(pos_ion_annotations)):
            if abs(pos_matrix[i][j]) > z_score_threshold:
                pos_affected_metabolites[gene_annotations[j]].append((i, pos_matrix[i][j]))

    cur.execute('select gene_id, accession from genes where '
                'strain_id = ANY(SELECT strain_id from strain WHERE strain.strain = %s)', ('mg1655', ))

    accession_to_gene_id = dict()
    for record in cur:
        accession_to_gene_id[record['accession'].upper()] = record['gene_id']

    upload_pos_tuples = []
    sql = prepare_tuple_style_sql_query('metabolomics', 'mg1655', SPECIES_SCHEMA['metabolomics'])
    for gene in pos_affected_metabolites:

        # this ignores pseudogenes
        converted_name = unified_std.convert(gene, 'Escherichia coli', strain='MG1655')
        if converted_name == gene.upper():
            continue
        if converted_name.upper() not in accession_to_gene_id:
            continue

        for (m, zscore) in pos_affected_metabolites[gene]:
            (ion_mass, kegg, name) = pos_ion_annotations[m]
            upload_pos_tuples.append((accession_to_gene_id[converted_name.upper()], ion_mass_to_pk[ion_mass]))

    psycopg2.extras.execute_values(cur, sql, argslist=upload_pos_tuples, page_size=2000)

    upload_neg_tuples = []
    for gene in neg_affected_metabolites:
        converted_name = unified_std.convert(gene, 'Escherichia coli', strain='MG1655')
        if converted_name == gene.upper():
            continue
        if converted_name.upper() not in accession_to_gene_id:
            continue
        for (m, zscore) in neg_affected_metabolites[gene]:
            (ion_mass, kegg, name) = neg_ion_annotations[m]
            upload_neg_tuples.append((accession_to_gene_id[converted_name.upper()], ion_mass_to_pk[ion_mass]))

    psycopg2.extras.execute_values(cur, sql, argslist=upload_neg_tuples, page_size=2000)



def build_metabolite_table(cur):

    """

    Built separately since it depends on files that are generated during DB uploads.

    :return:
    """
    build_sauer_metabolite_table(cur, os.path.join(constants.INPUT_DIR, 'metabolomics'))


if __name__ == '__main__':

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
        main(cur=cursor)

    connect.commit()
    connect.close()
