import glob
import logging
import os.path
from collections import defaultdict

import Bio.SeqIO
import psycopg2
import psycopg2.extras

from resistome import constants
from resistome.utils import database_utils
from resistome.sql.sql_build_utils import SPECIES_SCHEMA, prepare_sql_query, prepare_tuple_style_sql_query, \
    insert_genome_sequences, insert_sequence_data, insert_strain, parse_uniprot, \
    insert_uniprot_data, \
    insert_regulon_db_features, \
    parse_essential_genes_file, \
    build_mutational_prediction_table, SNAP2_GENES, parse_string_file
from resistome.utils.regulondb_parser import parse_regulondb_distribution
from typing import List, Tuple, Union, Dict, Set

log = logging.getLogger(__name__)


def fetch_target_file(target_directory: str, name_filter: str) -> str:

    """

    Extracts files from the target directory matching the provided name filter. wildcards are automatically
    inserted before and after the filter.

    A ValueError is thrown if too many files are found, and an AssertionError if no files are found.

    :param target_directory:
    :param name_filter:
    :return:
    """

    files_found = []
    for filepath in glob.glob(os.path.join(target_directory, '*' + name_filter + '*')):
        files_found.append(filepath)

    if len(files_found) > 1:
        raise ValueError('Too many files found! %s' % files_found)

    try:
        return files_found[0]
    except IndexError:
        raise AssertionError('Did not find any files matching the requested filter (%s) in %s'
                             % (name_filter, target_directory))


def extract_synonyms_from_genbank(genbank_file: str, proteome: str) -> Dict[str, Set[str]]:

    """

    Extracts extra synonyms from Genbank files. Depending on the format, this function will:

    Map synonyms to gene locus tags (accessions) or formal genes. Locus tags are used if available.

    Synonyms are found in the following fields:
        notes
        product
        gene

    And in the Uniprot tab file for the K-12 proteome. The actual gene synonym field is extracted elsewhere.

    :param genbank_file:
    :param proteome:
    :return:
    """

    synonyms = defaultdict(set)
    for record in Bio.SeqIO.parse(open(genbank_file), format='genbank'):

        for feature in record.features:
            if 'locus_tag' not in feature.qualifiers and 'gene' not in feature.qualifiers:
                continue

            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers['locus_tag'][0]
            else:
                locus_tag = feature.qualifiers['gene'][0]

            if 'gene' in feature.qualifiers and 'ins' in feature.qualifiers['gene'][0]:
                # get proper insertion element name
                if 'note' in feature.qualifiers and 'ins' in feature.qualifiers['note'][0] and ';' not in \
                        feature.qualifiers['note'][0]:
                    synonyms[locus_tag].add(feature.qualifiers['note'][0])
            if 'note' in feature.qualifiers and ':' in feature.qualifiers['note'][0] and 'similar to' not in \
                    feature.qualifiers['note'][0]:
                # get JW, ECK, or B numbers
                notes = feature.qualifiers['note'][0].split(';')
                for token in notes:
                    if 'ECK' in token or 'JW' in token or 'b' in token and 'ins' not in locus_tag:
                        # fuck you too buddy.
                        synonyms[locus_tag].update(filter(lambda x: ('b' in x and str.isnumeric(x[x.find('b') + 1]))
                                                                    or 'ECK' in x
                                                                    or 'JW' in x, token.split(':')))
            if 'note' in feature.qualifiers and (
                    'similar to' in feature.qualifiers['note'][0] or 'synonym' in feature.qualifiers['note'][0]):
                # also fuck you
                # get implicitly mapped tags
                for token in feature.qualifiers['note'][0].split('; '):
                    if ('similar to' in token or 'synonym' in token) and 'b' in token and str.isnumeric(
                            token[token.find('b') + 1]):
                        bnumber = token.replace('similar to', '').replace('synonym', '').strip()
                        if '(' in bnumber:
                            bnumber = bnumber[bnumber.find('(') + 1:bnumber.find(')')]
                        synonyms[locus_tag].add(bnumber)
            if 'product' in feature.qualifiers:
                # extract protein names
                product_string = feature.qualifiers['product'][0].split(' ')
                gene_name = None
                if len(product_string[-1]) == 4:
                    gene_name = product_string[-1]
                elif len(product_string[0].replace(',', '')) == 4:
                    gene_name = product_string[0].replace(',', '')

                if gene_name is not None and (str.isupper(gene_name[0]) or str.isupper(gene_name[-1])) and gene_name != 'Type':
                    gene_name = gene_name.strip().replace('-', '')
                    if len(gene_name.strip()) > 0:
                        synonyms[locus_tag].add(gene_name)

            synonyms[locus_tag].update(feature.qualifiers.get('gene_synonym', []))

            fixed_syns = set()
            for x in synonyms[locus_tag]:
                if ';' in x:
                    fixed_syns.update(x.split(';'))
                elif '+' in x:
                    fixed_syns.update(x.split('+'))
                else:
                    fixed_syns.add(x)

            synonyms[locus_tag.upper()] = fixed_syns
            if feature.location.strand == 1:
                strand = '+'
            else:
                strand = '-'
            # also add exact location mapping; +1 to start to compensate for zero-indexing
            synonyms[(str(int(feature.location.start + 1)), str(int(feature.location.end)), strand)] = fixed_syns

    if proteome is not None:
        with open(proteome) as f:
            # map uniprot standard names
            header = {x: k for k, x in enumerate(f.readline().strip().split('\t'))}
            for line in f:
                tokens = line.strip().split('\t')
                syns_of_our_fathers = tokens[header['Gene names']].split(' ')
                for syn in syns_of_our_fathers[1:]:
                    synonyms[syns_of_our_fathers[0].upper()].update(syn.split('/'))

    output = dict()
    for x, y in synonyms.items():
        if len(y) > 0 and len(x) > 0:
            output[x] = y
            assert '12' not in y

    return output


def extract_xrefs_from_genbank(genbank_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:

    db_xref = dict()
    ecocyc = dict()
    for record in Bio.SeqIO.parse(open(genbank_file), format='genbank'):

        for feature in record.features:
            if 'locus_tag' not in feature.qualifiers:
                continue

            locus_tag = feature.qualifiers['locus_tag'][0]

            if 'db_xref' in feature.qualifiers:

                for ref in feature.qualifiers['db_xref']:
                    if 'UniProt'.upper() in ref.upper():
                        uniprot_id = ref.split(':')[-1]
                        db_xref[uniprot_id] = locus_tag.upper()
                    elif 'ECOCYC' in ref.upper():
                        ecocyc_id = ref.split(':')[-1]
                        ecocyc[ecocyc_id] = locus_tag.upper()

    return db_xref, ecocyc


def parse_ncbi_feature_file(feature_file: str, genbank_file: str,
                            proteome_file: str,
                            additional_mappings: Dict[str, Set[str]] = None) -> List[Tuple[Tuple[Union[str, int, List], ...], Tuple[Union[str, int, List], ...]]]:

    """

    Handles parsing of NCBI feature files into the gene table expected in the Resistome support database. The
    Genbank file is used for synonym extraction.

    :param feature_file:
    :return:
    """

    if additional_mappings is None:
        additional_mappings = dict()

    seen_before = set()

    output_tuples = []

    essential_set = parse_essential_genes_file(os.path.join(constants.INPUT_DIR,
                                                            'biological_info',
                                                            'essential_genes_ecoli.txt'))

    synonyms = extract_synonyms_from_genbank(genbank_file=genbank_file,
                                             proteome=proteome_file)

    with open(feature_file) as f:
        header = {x: k for k, x in enumerate(f.readline().strip().split('\t'))}
        for line in f:
            tokens = line.strip().split('\t')
            entity_class = tokens[header['class']]

            if entity_class == 'protein_coding':
                continue

            locus_tag = tokens[header['locus_tag']].upper()
            start = tokens[header['start']]
            end = tokens[header['end']]
            direction = tokens[header['strand']]
            protein_accession = tokens[header['product_accession']]
            description = tokens[header['name']]

            name = tokens[header['symbol']]

            if (locus_tag, name) in seen_before:
                continue
            seen_before.add((locus_tag, name))

            if name == '':
                name = locus_tag

            if description == '':
                description = None

            assert len(locus_tag) > 0

            # ['biocyc_id',
            #  'name',
            #  'synonyms',
            #  'products',
            #  'accession',
            #  'start',
            #  'stop',
            #  'direction',
            #  'essential']

            synonym_upload = []
            if locus_tag.upper() in synonyms:
                synonym_upload.extend(list(synonyms[locus_tag.upper()]))
            if name.upper() in synonyms:
                synonym_upload.extend(list(synonyms[name.upper()]))
            if locus_tag.upper() in additional_mappings:
                synonym_upload.extend(additional_mappings[locus_tag.upper()])
            if name.upper() in additional_mappings:
                synonym_upload.extend(additional_mappings[name.upper()])
            if (start, end, direction) in additional_mappings:
                synonym_upload.extend(additional_mappings[(start, end, direction)])
            synonym_upload = {x.upper() for x in set(synonym_upload)}

            synonym_upload = list(filter(lambda x: len(x) > 0, synonym_upload))

            gene_tuple = (locus_tag.upper(), name.upper(), synonym_upload, [protein_accession], locus_tag.upper(),
                          description,
                          'pseudogene' in entity_class,
                          name.upper() in essential_set)

            position_tuple = (int(start), int(end), direction)
            output_tuples.append((gene_tuple, position_tuple))

    return output_tuples


def get_go_terms(gene_accession_tuples: List[Tuple[str, str]], go_file = None):

    if go_file is None:
        go_file = os.path.join(constants.INPUT_DIR, 'ontologies', 'ecocyc_DEC2020.gaf')

    gene_to_GO = defaultdict(set)
    output_tuples = []

    with open(go_file) as f:
        for line in f:
            if line[0] == '!':
                # comment
                continue

            tokens = line.strip().split('\t')
            gene = tokens[2].upper()
            go_term = tokens[4]
            gene_to_GO[gene].add(go_term)

    for (accession, gene) in gene_accession_tuples:
        if gene.upper() in gene_to_GO:
            output_tuples.extend([(accession, x) for x in gene_to_GO[gene.upper()]])

    return output_tuples


def main(cur: psycopg2._psycopg.cursor, source_data: str, species: str, strain: str):

    """

    Handles the main parsing for NCBI and NCBI adjacent files.

    :param cur:
    :param source_data:
    :param species:
    :param strain:
    :return:
    """

    strain_id = insert_strain(cur, species, strain)

    # TODO parse feature table (gff?)
    # TODO cross reference sequence files
    feature_file = fetch_target_file(source_data, 'feature_table')

    try:
        genome = fetch_target_file(source_data, 'v1_genomic.fna')
    except AssertionError:
        # v2 is sometimes used by NCBI
        genome = fetch_target_file(source_data, 'v2_genomic.fna')

    cds_sequences = fetch_target_file(source_data, 'cds_from_genomic')
    protein_sequences = fetch_target_file(source_data, 'protein')
    rna_sequences = fetch_target_file(source_data, 'rna_from_genomic')

    if strain != 'w3110':
        genbank = fetch_target_file(source_data, 'gbff')
        gca_mapping = None
    else:
        # RefSeq has the accession, Genbank has the mappings...
        # why?????????????????????????????????????????????????????????????????
        # technology was a mistake.
        # need to special case this issue.
        genbank = fetch_target_file(source_data, 'GCF_000010245.2_ASM1024v1_genomic.gbff')
        gca_mapping = fetch_target_file(source_data, 'GCA_000010245.1_ASM1024v1_genomic.gbff')

    # this is the K-12 protein mapping used to identify extra synonyms for E. coli genes
    proteome_mapping = fetch_target_file(os.path.join(constants.INPUT_DIR, 'biocyc', 'generic'),
                                         'UP000000625.tab')

    # get xrefs to enable indexing other database data
    uniprot_cross_refs, ecocyc_to_accession = extract_xrefs_from_genbank(genbank)

    if len(uniprot_cross_refs.keys()) > 0:
        # if you have any uniprot refs, we can use this data from uniprot to associate your genes
        # with uniprot protein annotations.
        # granted, we could probably get gffs for ech organism instead...
        uniprot_to_annotations = parse_uniprot(os.path.join(constants.INPUT_DIR, 'biocyc', 'generic',
                                                            'uniprot-proteome_UP000000625.gff'))
    else:
        uniprot_to_annotations = dict()

    if gca_mapping is not None:
        # pull out extra IDs for W3110
        result = extract_synonyms_from_genbank(gca_mapping, proteome_mapping)
    else:
        # all other cases
        result = dict()

    # parses the feature table to get positions for all genes
    # these tuples are directly uploaded into the db later
    gene_features = parse_ncbi_feature_file(feature_file, genbank, proteome_mapping,
                                            additional_mappings=result)

    genome_dict = dict()
    # this is fine for multi part genomes.
    for record in Bio.SeqIO.parse(open(genome), format='fasta'):
        genome_dict[record.id] = str(record.seq)

    protein_dict = dict()
    cds_dict = dict()
    protein_id_to_locus = dict()
    for record in Bio.SeqIO.parse(open(cds_sequences), format='fasta'):

        description = record.description
        tokens = description.split(' ')
        target_token = None
        for token in tokens:
            if 'locus_tag' in token:
                target_token = token
        # assert every CDS has a unique acccesion
        assert target_token is not None and 'locus_tag' in target_token, target_token
        locus_tag = target_token.replace('[', '').replace(']', '').split('=')[1].upper()

        nt_seq = str(record.seq)
        cds_dict[locus_tag] = nt_seq
        if 'pseudo=true' not in description:
            # not a pseudogene, get protein ID
            target_token = None
            for token in tokens:
                if 'protein_id' in token:
                    target_token = token

            protein_id = target_token.replace('[', '').replace(']', '').split('=')[1]
            protein_id_to_locus[protein_id] = locus_tag.upper()

    for record in Bio.SeqIO.parse(open(rna_sequences), format='fasta'):

        # same thing but for genomic RNA
        description = record.description
        tokens = description.split(' ')
        target_token = None
        for token in tokens:
            if 'locus_tag' in token:
                target_token = token
        assert target_token is not None and 'locus_tag' in target_token
        locus_tag = target_token.replace('[', '').replace(']', '').split('=')[1]
        nt_seq = str(record.seq)
        cds_dict[locus_tag.upper()] = nt_seq

    for record in Bio.SeqIO.parse(open(protein_sequences), format='fasta'):
        # prot sequences.
        # note that the protein ids are used later to associate DeMaSk output to the original gene locus.
        locus_tag = protein_id_to_locus[record.id]
        protein_dict[locus_tag.upper()] = str(record.seq)

    # get existing GO ids
    cur.execute('select go_id, go_term from go_table')

    go_term_to_go_id = dict()
    for record in cur:
        go_term_to_go_id[record['go_term']] = record['go_id']

    # should do this in batch but I don't really care...
    accession_to_gene_id = dict()
    SQL = prepare_tuple_style_sql_query('genes', strain, SPECIES_SCHEMA['genes']) + ' RETURNING gene_id, accession'
    gene_tuples = []
    for (gene_feature, pos_feature) in gene_features:
        # TODO optimize with batch upload
        # uploads all genes, locations
        gene_tuples.append((strain_id,) + gene_feature)

    results = psycopg2.extras.execute_values(cur, SQL, argslist=gene_tuples, page_size=4000, fetch=True)
    for result in results:
        accession_to_gene_id[result['accession'].upper()] = result['gene_id']

    pos_tuples = []
    for (gene_feature, pos_feature) in gene_features:

        gene_id = accession_to_gene_id[gene_feature[0].upper()]
        pos_tuples.append((gene_id,) + pos_feature)

    psycopg2.extras.execute_values(cur,
                                   sql='INSERT INTO gene_locations (gene_id, start, stop, direction) VALUES %s',
                                   argslist=pos_tuples, page_size=4000)

    insert_genome_sequences(cur=cur, schema=strain, columns=SPECIES_SCHEMA['genome'],
                            strain_id=strain_id,
                            genome_object=genome_dict)

    nt_tuples = []
    sql = prepare_tuple_style_sql_query('nt_sequences', strain, SPECIES_SCHEMA['dna'])
    for locus_tag, seq in cds_dict.items():
        nt_tuples.append((accession_to_gene_id[locus_tag.upper()], seq))

    psycopg2.extras.execute_values(cur, sql, nt_tuples, page_size=2000)

    aa_tuples = []
    sql = prepare_tuple_style_sql_query('aa_sequences', strain, SPECIES_SCHEMA['aa'])
    for locus_tag, seq in protein_dict.items():
        aa_tuples.append((accession_to_gene_id[locus_tag.upper()], seq))

    psycopg2.extras.execute_values(cur, sql, aa_tuples, page_size=2000)

    if len(uniprot_to_annotations.keys()) > 0:
        insert_uniprot_data(cur=cur, schema=strain, columns=SPECIES_SCHEMA['uniprot'],
                            uniprot_dict=uniprot_to_annotations, uniprot_to_accession=uniprot_cross_refs,
                            accession_to_gene_id=accession_to_gene_id)

    if strain == 'mg1655':
        # hard-coding this since I don't know another way to do it...
        # insert the genomic annotations from RegulonDB, including the regulatory networks.
        # these are parsed in the imported functions from the regulondb parser. Source files are included in the
        # repo. Current version: 10.6
        insert_regulon_db_features(accession_to_gene_id=accession_to_gene_id,
                                   database_cursor=cur,
                                   schema=strain,
                                   strain_id=strain_id)

        regulatory_interactions= parse_regulondb_distribution(strain_id=strain_id,
                                                              accession_to_gene_id=accession_to_gene_id)

        sql = prepare_tuple_style_sql_query('interactions', strain, SPECIES_SCHEMA['interactions'])
        psycopg2.extras.execute_values(cur, sql, argslist=regulatory_interactions, page_size=2000)

        strings_db = os.path.join(constants.INPUT_DIR, 'interactions', '511145.protein.links.v11.0.txt.gz')
        # the docs do not provide any guidance on select this beyond a single example where score = 0.4: medium
        # confidence. I selected this threshold arbitrarily, so it might be too restrictive.
        protein_edges = parse_string_file(strings_db, minimum_score=0.70)
        upload_tuples = []
        for (p1, p2) in protein_edges:
            upload_tuples.append((strain_id, 'protein-protein',
                                  p1, p2,
                                  accession_to_gene_id[p1], accession_to_gene_id[p2],
                                  '?',
                                  'strings'))

        psycopg2.extras.execute_values(cur, sql=sql, argslist=upload_tuples, page_size=5000)

    # add in predictions from these tools
    # You can easily change this function to load data for every gene, but it is a titanic amount of data.
    build_mutational_prediction_table(cur=cur, strain_to_process=strain,
                                      unique_id_to_accession_dict=protein_id_to_locus,
                                      accession_to_gene_id=accession_to_gene_id,
                                      methods={'inps', 'snap2', 'demask'},
                                      variant_effect_predictor_genes=SNAP2_GENES)

    go_terms = get_go_terms([(x[0][0], x[0][1]) for x in gene_features])
    go_tuples = []
    SQL = prepare_tuple_style_sql_query('go_terms', strain_id, SPECIES_SCHEMA['go_terms'])

    for locus_tag, go_term in go_terms:
        go_tuples.append((accession_to_gene_id[locus_tag], go_term_to_go_id[go_term]))

    psycopg2.extras.execute_values(cur, SQL, argslist=go_tuples, page_size=2000)


if __name__ == '__main__':
    # todo: add table of dna binding proteins - consensus sequences

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    root_dir = constants.INPUT_DIR
    ECOLI_BW25113_DIR = os.path.join(root_dir, 'biocyc', 'bw25113')
    ECOLI_W3110_DIR = os.path.join(root_dir, 'biocyc', 'w3110')
    ECOLI_MDS42_DIR = os.path.join(root_dir, 'biocyc', 'mds42')
    ECOLI_BL21_DIR = os.path.join(root_dir, 'biocyc', 'bl21')
    ECOLI_BL21_DE3_DIR = os.path.join(root_dir, 'biocyc', 'bl21_de3_')

    ECOLI_MG1655_DIR = os.path.join(root_dir, 'biocyc', 'mg1655', 'data')
    ECOLI_B_DIR = os.path.join(root_dir, 'biocyc', 'rel606', 'data')
    ECOLI_W_DIR = os.path.join(root_dir, 'biocyc', 'ecoliw', 'data')

    # source_dirs = [ECOLI_MG1655_DIR, ECOLI_W3110_DIR, ECOLI_BW25113_DIR, ECOLI_MDS42_DIR, ECOLI_BL21_DIR,
    #                ECOLI_BL21_DE3_DIR,
    #                ECOLI_B_DIR, ECOLI_W_DIR]

    source_dirs = [ECOLI_B_DIR]

    ncbi_strain_tuples = [(x, os.path.split(x)[1]) for x in source_dirs]

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:

        with open(os.path.join(constants.INPUT_DIR, 'sql', 'resistome_biocyc_schema.sql'), 'r') as f:
            sql_schema = ''.join(f.readlines())
            # make schema species specific
            cursor.execute(sql_schema)

        for source_data, strain in ncbi_strain_tuples:
            main(cursor, source_data, 'Escherichia coli', strain)

    connect.commit()
    connect.close()
