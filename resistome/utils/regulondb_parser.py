from resistome.constants import INPUT_DIR
import os
from typing import List, Tuple, Optional, Dict
from collections import defaultdict
import itertools
import json


REGULON_DB_DIR = os.path.join(INPUT_DIR, 'regulondb')


def load_tu_groups(tu_filepath: str) -> Dict[str, List[str]]:

    if not os.path.exists(tu_filepath):
        raise AssertionError('Did not find %s' % tu_filepath)

    tu_to_genes = defaultdict(list)
    with open(tu_filepath) as f:
        for line in f:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            tu_id = tokens[0]
            gene_id = tokens[1]

            tu_to_genes[tu_id].append(gene_id)

    return  tu_to_genes


def load_operon_groups(operon_filepath: str, gene_starts: Dict[int, str], gene_ends: Dict[int, str]) -> Dict[str, List[str]]:

    """

    # 1) OBJECT_ID
    # 2) OBJECT_NAME
    # 3) OBJECT_TYPE
    # 4) OBJECT_STRAND
    # 5) OBJECT_POSLEFT
    # 6) OBJECT_POSRIGHT
    # 7) OBJECT_COLOR
    # 8) TOOL_TIP
    # 9) LINE_TYPE
    # 10) LABEL_SIZE

    :param operon_filepath:
    :param gene_starts:
    :param gene_ends:
    :return:
    """

    if not os.path.exists(operon_filepath):
        raise AssertionError('Did not find %s' % operon_filepath)

    operon_to_genes = dict()
    with open(operon_filepath) as f:
        for line in f:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            object_type = tokens[2]

            if object_type == 'operon':
                operon_accession = tokens[0]
                operon_start = int(tokens[4])
                operon_end = int(tokens[5])

                # get all genes in between the operon start and end
                # these are the genes associated with the operon
                genes_included = [gene_starts[x] for x in filter(lambda y: operon_start <= y <= operon_end,
                                                                 gene_starts.keys())]

                operon_to_genes[operon_accession] = genes_included

    return operon_to_genes


def load_gene_accessions(gene_filepath: str) -> Tuple[Dict[int, str], Dict[int, str]]:

    # Columns:
    # 1) GENE_ID
    # 2) GENE_NAME
    # 3) GENE_POSLEFT
    # 4) GENE_POSRIGHT
    # 5) GENE_STRAND
    # 6) GENE_SEQUENCE
    # 7) GC_CONTENT
    # 8) CRI_SCORE
    # 9) GENE_NOTE
    # 10) GENE_INTERNAL_COMMENT
    # 11) KEY_ID_ORG
    # 12) GENE_TYPE

    if not os.path.exists(gene_filepath):
        raise AssertionError('Did not find %s' % gene_filepath)

    starts = dict()
    ends = dict()

    with open(gene_filepath) as f:
        for line in f:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            accession = tokens[0]

            if tokens[2] == '':
                continue

            start = int(tokens[2])
            end = int(tokens[3])

            starts[start] = accession
            ends[end] = accession

    return starts, ends


def load_regulondb_external_ids(external_id_filepath: str, target_col: int = 2) -> Dict[str, str]:

    """

    These are internal accession => b-number mapping.

    :param external_id_filepath:
    :return:
    """

    if not os.path.exists(external_id_filepath):
        raise AssertionError('Did not find %s' % external_id_filepath)

    output_dict = dict()
    with open(external_id_filepath) as f:
        for line in f:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            # if it is a K-12 accession and the external id is a b-number
            if 'ECK12' in tokens[0] and tokens[target_col][0].upper() == 'B':
                rgdb_accession = tokens[0]
                bnumber = tokens[target_col]
                output_dict[rgdb_accession] = bnumber.upper()

    return output_dict


def load_generic_mapping(sigma_filepath, gene_position: int = 3) -> Dict[str, List[str]]:

    """

    # # RegulonDB: 10.6
    # Table: SIGMA_TMP
    # Columns:
    # 1) SIGMA_ID
    # 2) SIGMA_NAME
    # 3) SIGMA_SYNONYMS
    # 4) SIGMA_GENE_ID
    # 5) SIGMA_GENE_NAME
    # 6) SIGMA_COREGULATORS
    # 7) SIGMA_NOTES
    # 8) SIGMA_SIGMULON_GENES
    # 9) KEY_ID_ORG

    :param sigma_filepath:
    :return:
    """

    if not os.path.exists(sigma_filepath):
        raise AssertionError('Did not find %s' % sigma_filepath)

    mapping = defaultdict(list)
    with open(sigma_filepath) as f:
        for line in f:
            if line[0] == '!' or line[0] == '#':
                continue

            tokens = line.strip().split('\t')
            # print(tokens)
            sigma_id = tokens[0]
            gene_id = tokens[gene_position]

            mapping[sigma_id].append(gene_id)

    return mapping


def parse_regulondb_distribution(strain_id, accession_to_gene_id) -> List[Tuple[str, ...]]:

    def id_disambiguator(accession: str) -> List[str]:

        if accession in eck_to_bnumber:
            return [accession]
        if accession in tf_mappings:
            return tf_mappings[accession]
        if accession in tu_to_gene:
            return tu_to_gene[accession]
        if accession in operon_to_gene:
            return operon_to_gene[accession]
        if accession in product_mappings:
            return product_mappings[accession]
        if accession in sigma_mappings:
            return sigma_mappings[accession]

        return []

    tu_file = os.path.join(REGULON_DB_DIR, 'tu_gene_link.txt')
    tf_file = os.path.join(REGULON_DB_DIR, 'product_tf_link.txt')
    sigma_file = os.path.join(REGULON_DB_DIR, 'sigma_tmp.txt')
    product_file = os.path.join(REGULON_DB_DIR, 'gene_product_link.txt')
    product_info = os.path.join(REGULON_DB_DIR, 'product.txt')
    object_props = os.path.join(REGULON_DB_DIR, 'objects_properties_tmp.txt')
    external_refs = os.path.join(REGULON_DB_DIR, 'object_external_db_link.txt')
    synonyms = os.path.join(REGULON_DB_DIR, 'object_synonym.txt')
    genes = os.path.join(REGULON_DB_DIR, 'gene.txt')

    tu_to_gene = load_tu_groups(tu_file)
    gene_starts, gene_ends = load_gene_accessions(genes)
    sigma_mappings = load_generic_mapping(sigma_file)
    tf_mappings = load_generic_mapping(tf_file, gene_position=1)
    product_mappings = load_generic_mapping(product_file, gene_position=1)
    product_info_last_ditch = load_generic_mapping(product_info, gene_position=1)
    gene_info_last_dict = load_generic_mapping(genes, gene_position=1)
    operon_to_gene = load_operon_groups(object_props, gene_starts, gene_ends)
    regulatory_network = os.path.join(REGULON_DB_DIR, 'genetic_network.txt')

    eck_to_bnumber = load_regulondb_external_ids(external_refs)
    eck_to_synonym = load_regulondb_external_ids(synonyms, target_col=1)

    for eck in eck_to_synonym:
        if eck not in eck_to_bnumber:
            eck_to_bnumber[eck] = eck_to_synonym[eck]

    # 1) REGULATOR_ID
    # 2) REGULATOR_NAME
    # 3) REGULATED_ID
    # 4) REGULATED_NAME
    # 5) FUNCTION_INTERACTION
    # 6) EVIDENCE
    # 7) REGULATOR_TYPE
    # 8) REGULATED_TYPE

    regulator_type = 6
    regulated_type = 7
    interaction_type = 4

    a_id = 0
    b_id = 2

    output_interactions = []
    with open(regulatory_network) as f:
        for line in f:
            if line[0] == '#':
                continue

            tokens = line.strip().split('\t')

            regulator_accession = tokens[a_id]
            regulated_accession = tokens[b_id]
            interaction_dir = tokens[interaction_type]

            a_type = tokens[regulator_type]
            b_type = tokens[regulated_type]

            flattened_regulators = id_disambiguator(regulator_accession)
            flattened_regulated = id_disambiguator(regulated_accession)

            final_regulators = []
            final_regulated = []

            for regulator in flattened_regulators:
                if regulator in eck_to_bnumber:
                    final_regulators.append(eck_to_bnumber[regulator])
                elif regulator in product_info_last_ditch:
                    final_regulators.extend(product_info_last_ditch[regulator])
                elif regulator in gene_info_last_dict:
                    final_regulators.extend(gene_info_last_dict[regulator])
                else:
                    raise AssertionError

            if len(final_regulators) == 0:
                continue

            for regulated in flattened_regulated:

                if regulated in eck_to_bnumber:
                    final_id = [eck_to_bnumber[regulated]]
                elif regulated in product_info_last_ditch:
                    final_id = product_info_last_ditch[regulated]
                elif regulated in gene_info_last_dict:
                    final_id = gene_info_last_dict[regulated]
                else:
                    raise AssertionError

                final_regulated.extend(final_id)

            if len(final_regulated) == 0:
                continue

            actual_regulator = final_regulators[0]

            interaction_name = '%s-%s' % (a_type, b_type)
            for regulator, regulated in itertools.product(final_regulators, final_regulated):
                output_interactions.append((strain_id,
                                            interaction_name,
                                            actual_regulator.upper(),
                                            regulated.upper(),
                                            accession_to_gene_id.get(actual_regulator.upper(), None),
                                            accession_to_gene_id.get(regulated.upper(), None),
                                            interaction_dir,
                                            'regulondb'))

    return output_interactions


def generic_file_parser(filepath: str,
                        positions_to_tuple: List[str],
                        positions_to_dict: List[str],
                        tuple_placeholder: Tuple[int, str] = None,
                        exclude=None) -> List[Tuple[List[str], Dict[str, str]]]:

    if not os.path.exists(filepath):
        raise AssertionError('Did not find %s' % filepath)

    col_to_pos = dict()
    pos_to_col = dict()
    with open(filepath) as f:
        for line in f:
            if line[0] == '!' or line[0] == '#':
                if 'Columns:' in line:
                    line = f.readline()
                    while line[0] == '#' or line[0] == '!':
                        tokens = line.strip().split(' ')
                        position = int(tokens[1].replace(')', '')) - 1
                        name = tokens[-1]
                        col_to_pos[name] = position
                        pos_to_col[position] = name
                        line = f.readline()
                    break
            else:
                break

    output_tuples = []
    with open(filepath) as f:
        for line in f:
            if line[0] == '!' or line[0] == '#' or line[0] == '<':
                continue

            row_dict = dict()
            tokens = line.strip().split('\t')
            row = []

            skip = False
            if exclude is not None:
                for column, bool_function in exclude.items():
                    if not bool_function(tokens[col_to_pos[column]]):
                        skip = True
                        break
            if skip:
                continue

            for k in positions_to_tuple:
                pos = col_to_pos[k]
                token = tokens[pos]
                row.append(token)

            if tuple_placeholder is not None:
                row.insert(tuple_placeholder[0], tuple_placeholder[1])

            for k in positions_to_dict:
                pos = col_to_pos[k]
                row_dict[k] = tokens[pos]

            output_tuples.append((tuple(row), row_dict))

    return output_tuples


def handle_dna_binding_site_association(eck_to_bnumber: Dict[str, str]):

    dna_binding_sites = os.path.join(REGULON_DB_DIR, 'site.txt')
    regulatory_network = os.path.join(REGULON_DB_DIR, 'genetic_network.txt')
    reg_interactions = os.path.join(REGULON_DB_DIR, 'regulatory_interaction.txt')
    tf_gene_interactions = os.path.join(REGULON_DB_DIR, 'tf_gene_interaction.txt')
    conformation = os.path.join(REGULON_DB_DIR, 'conformation.txt')
    tf_file = os.path.join(REGULON_DB_DIR, 'product_tf_link.txt')
    tf_mappings = load_generic_mapping(tf_file, gene_position=1)
    mapping_to_tf = dict()
    for x, y in tf_mappings.items():
        for z in y:
            mapping_to_tf[z] = x

    gene_info_last_dict = load_generic_mapping(regulatory_network, gene_position=1)
    gene_to_eck = dict()

    for x, y in gene_info_last_dict.items():
        for z in y:
            gene_to_eck[z] = x

    site_to_regulator = defaultdict(list)
    with open(reg_interactions) as f:
        for line in f:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            site_id = tokens[3]
            regulator_maybe_id = tokens[1]
            site_to_regulator[site_id].append(regulator_maybe_id)

    with open(tf_gene_interactions) as f:
        for line in f:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            site_id = tokens[3]
            regulator_maybe_id = tokens[1]
            site_to_regulator[site_id].append(regulator_maybe_id)

    conformation_to_tf = defaultdict(list)
    with open(conformation) as f:
        for line in f:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            conformation_id = tokens[0]
            regulator_maybe_id = tokens[1]
            conformation_to_tf[conformation_id] = regulator_maybe_id

    binding_sites = generic_file_parser(dna_binding_sites,
                                        ['SITE_ID', 'SITE_POSLEFT', 'SITE_POSRIGHT'],
                                        ['SITE_SEQUENCE'],
                                        tuple_placeholder=(1, 'dna_binding_site'))

    temp_bs = []
    for (row, row_dict) in binding_sites:

        site_id = row[0]

        if site_id in site_to_regulator:
            bnumbers = []
            for conformation_id in site_to_regulator[site_id]:
                tf_id = conformation_to_tf[conformation_id]

                if tf_id in eck_to_bnumber:
                    bnumbers.append(eck_to_bnumber[tf_id])
                else:
                    for tf_complex_component_id in tf_mappings[tf_id]:
                        bnumbers.append(eck_to_bnumber[tf_complex_component_id])

            row_dict['tf_binding'.upper()] = bnumbers
        else:
            row_dict['tf_binding'.upper()] = []
        temp_bs.append((row, row_dict))

    return temp_bs


def handle_promoter_extraction(promoter_path):

    promoters_sites = generic_file_parser(promoter_path,
                                          ['PROMOTER_ID', 'POS_1'],
                                          ['PROMOTER_NAME', 'SIGMA_FACTOR', 'PROMOTER_SEQUENCE', 'PROMOTER_STRAND'],
                                          tuple_placeholder=(1, 'promoter'))

    sigma_to_bnumber = {'Sigma70': 'B3067', 'Sigma38': 'B2741',
                       'Sigma54': 'B3202', 'Sigma32': 'B3461',
                       'Sigma24': 'B2573', 'Sigma28': 'B1922'}

    output = []
    for (row, row_dict) in promoters_sites:

        if row[2] == '':
            continue

        position = int(row[2])
        rpos = position + len(row_dict['PROMOTER_SEQUENCE'])
        row = row + (rpos,)

        sigma = row_dict['SIGMA_FACTOR']

        if sigma != 'unknown' and len(sigma) > 0:
            sigma_tokens = [x.strip() for x in sigma.split(',')]
            row_dict['SIGMA_FACTOR'] = [sigma_to_bnumber[x] for x in sigma_tokens]
        else:
            row_dict['SIGMA_FACTOR'] = []

        output.append((row, row_dict))

    return output


def prepare_final_output(entry_json_tuples: List[Tuple[List[str], Dict[str, str]]], strain_id: int):

    output = []

    duplicate_id = set()
    for (row, row_dict) in entry_json_tuples:

        if row[0] in duplicate_id:
            duplicate_id.add(row[0])
            continue

        duplicate_id.add(row[0])

        new_annotation_dict = dict()
        for x, y in row_dict.items():
            if 'SEQUENCE' in x:
                new_annotation_dict['SEQUENCE'] = y
            elif 'STRAND' in x:
                new_annotation_dict['STRAND'] = '+' if y == 'forward' else '-'
            else:
                new_annotation_dict[x] = y

        row = (strain_id, ) + row + (json.dumps(new_annotation_dict), 'regulondb')

        try:
            assert len(row) == 7, row
        except AssertionError:
            print('Row did not have the expected number of columns.')
            raise
        output.append(row)

    return output


def extract_mg1655_genome_features(strain_id: int):

    """

    Loads genome annotations contained within RegulonDB; this is mainly meant to help supplement (replace?)
    Biocyc data which is difficult to parse and hard to access.

    :return:
    """

    # anti-terminators
    anti_term = os.path.join(REGULON_DB_DIR, 'attenuator_terminator.txt')
    # terminators
    terminator = os.path.join(REGULON_DB_DIR, 'terminator.txt')
    # riboswitches
    riboswitch = os.path.join(REGULON_DB_DIR, 'rfam.txt')
    # dna binding sites (I think?)
    # relegated to the super complex method above
    # SD boxes
    sd_boxes = os.path.join(REGULON_DB_DIR, 'shine_dalgarno.txt')
    # promoters
    promoters = os.path.join(REGULON_DB_DIR, 'promoter.txt')
    # sRNAs
    srnas_file = os.path.join(REGULON_DB_DIR, 'srna_interaction.txt')
    # operons, TSSs, etc
    object_props = os.path.join(REGULON_DB_DIR, 'objects_properties_tmp.txt')

    # help convert to standard IDs (for the sigma factors)
    external_refs = os.path.join(REGULON_DB_DIR, 'object_external_db_link.txt')
    synonyms = os.path.join(REGULON_DB_DIR, 'object_synonym.txt')

    eck_to_bnumber = load_regulondb_external_ids(external_refs, target_col=2)
    eck_to_synonym = load_regulondb_external_ids(synonyms, target_col=1)
    eck_to_name = dict()

    for eck in eck_to_synonym:
        if eck not in eck_to_bnumber:
            eck_to_bnumber[eck] = eck_to_synonym[eck]

    genes = os.path.join(REGULON_DB_DIR, 'gene.txt')
    gene_starts, gene_ends = load_gene_accessions(genes)

    with open(genes) as f:
        for line in f:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            eck_to_name[tokens[0]] = tokens[1]

    # create table <species>.genomic_features (
    #
    # feature_id serial primary key,
    # biocyc_id	text unique not null,
    # feature_type	text not null,
    # start	integer not null,
    # stop	integer not null,
    # associated_data jsonb,
    # source  text not null
    #
    # );

    anti_terminator_tuples = generic_file_parser(anti_term,
                                                 ['A_TERMINATOR_ID', 'A_TERMINATOR_TYPE', 'A_TERMINATOR_POSLEFT', 'A_TERMINATOR_POSRIGHT'],
                                                 ['A_TERMINATOR_SEQUENCE', 'A_TERMINATOR_ENERGY'])
    terminator_tuples = generic_file_parser(terminator,
                                            ['TERMINATOR_ID', 'TERMINATOR_POSLEFT', 'TERMINATOR_POSRIGHT'],
                                            ['TERMINATOR_CLASS', 'TERMINATOR_SEQUENCE'],
                                            tuple_placeholder=(1, 'terminator'))
    riboswitch_tuples = generic_file_parser(riboswitch,
                                     ['RFAM_ID', 'RFAM_POSLEFT', 'RFAM_POSRIGHT'],
                                     ['RFAM_TYPE', 'RFAM_SEQUENCE', 'RFAM_STRAND', 'RFAM_SCORE'],
                                     tuple_placeholder=(1, 'riboswitch'))

    dna_binding_sites = handle_dna_binding_site_association(eck_to_bnumber)

    sd_box_objects = generic_file_parser(sd_boxes,
                                         ['SHINE_DALGARNO_ID', 'SHINE_DALGARNO_POSLEFT', 'SHINE_DALGARNO_POSRIGHT'],
                                         ['SHINE_DALGARNO_DIST_GENE', 'SHINE_DALGARNO_SEQUENCE'],
                                         tuple_placeholder=(1, 'shine_dalgarno'))

    promoter_sites = handle_promoter_extraction(promoters)

    tss_entries = generic_file_parser(object_props,
                                      ['OBJECT_ID', 'OBJECT_POSLEFT', 'OBJECT_POSRIGHT'],
                                      ['OBJECT_STRAND', 'OBJECT_NAME'],
                                      tuple_placeholder=(1, 'tss'),
                                      exclude={'OBJECT_NAME': lambda x: 'TSS_' in x})

    operon_entries = generic_file_parser(object_props,
                                         ['OBJECT_ID', 'OBJECT_POSLEFT', 'OBJECT_POSRIGHT'],
                                         ['OBJECT_STRAND', 'OBJECT_NAME'],
                                         tuple_placeholder=(1, 'operon'),
                                         exclude={'OBJECT_TYPE': lambda x: x == 'operon'})

    for row, row_dict in operon_entries:

        start = int(row[-2])
        end = int(row[-1])

        containe_gene_starts = filter(lambda x: start <= x <= end, gene_starts)
        contained_genes = []
        for gene_start in containe_gene_starts:
            contained_genes.append(gene_starts[gene_start])

        operon_definition = []
        for eck in contained_genes:
            if eck in eck_to_bnumber:
                operon_definition.append(eck_to_bnumber[eck])
            else:
                operon_definition.append(eck_to_name[eck])

        row_dict['OPERON_GENES'] = operon_definition

    srnas = generic_file_parser(object_props,
                                ['OBJECT_ID', 'OBJECT_POSLEFT', 'OBJECT_POSRIGHT'],
                                ['OBJECT_STRAND', 'OBJECT_NAME'],
                                tuple_placeholder=(1, 'srnas'),
                                exclude={'OBJECT_TYPE': lambda x: x == 'srna'})

    srna_to_seq_target = dict()
    with open(srnas_file) as f:
        for line in f:
            if line[0] == '#':
                continue
            tokens = line.strip().split('\t')
            srna_to_seq_target[tokens[0].strip()] = (tokens[2], tokens[4])

    for row, row_dict in srnas:

        srna_id = row[0]
        target_eck, direction = srna_to_seq_target[srna_id]

        row_dict['TARGET'] = [eck_to_bnumber[target_eck]]

        if direction == 'activator':
            row_dict['EFFECT'] = '+'
        elif direction == 'repressor':
            row_dict['EFFECT'] = '-'
        elif direction == 'unknown' or len(direction) == 0:
            row_dict['EFFECT'] = "?"
        else:
            raise AssertionError('Direction %s' % direction)

    combined_output = []
    combined_output.extend(srnas)
    combined_output.extend(operon_entries)
    combined_output.extend(tss_entries)
    combined_output.extend(promoter_sites)
    combined_output.extend(sd_box_objects)
    combined_output.extend(dna_binding_sites)
    combined_output.extend(riboswitch_tuples)
    combined_output.extend(terminator_tuples)
    combined_output.extend(anti_terminator_tuples)

    insert_ready_output = prepare_final_output(combined_output, strain_id)

    print('Extracted %i genomic annotations from RegulonDB bundle' % len(insert_ready_output))

    gene_relationships = {'operon': 'OPERON_GENES',
                          'srnas': 'TARGET',
                          'promoter': 'SIGMA_FACTOR',
                          'dna_binding_site': 'TF_BINDING'}

    return insert_ready_output, gene_relationships


if __name__ == '__main__':
    extract_mg1655_genome_features()
