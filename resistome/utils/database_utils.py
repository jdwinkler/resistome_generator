import os
from collections import defaultdict
from resistome.sql.object_representation import ContainerClass
from resistome import constants
from resistome.utils.input_standardization import Standard


# expected input terms associated with paper, mutant, mutation, and transcriptome objects
def parse_input_terms():

    from collections import defaultdict

    input_terms = os.path.join(constants.INPUT_DIR,
                               'settings',
                               'InputFields_typed.txt')

    term_list_dict = defaultdict(list)
    variable_types = dict()
    parser_rules = dict()

    with open(input_terms, 'r') as f:

        headers = f.readline().strip().split('\t')
        c = {}

        for i, token in zip(range(0, len(headers)), headers):
            c[token] = i

        for line in f:
            line = line.translate(str.maketrans({'\n': None,
                                                 '\"': None}))
            tokens = line.split("\t")

            category = tokens[c['Category']]
            entry_name = tokens[c['Internal Name']]
            input_type = tokens[c['Type']]
            input_required = tokens[c['Required?']]
            multiple_values = tokens[c['Can select multiple?']]
            input_entry = tokens[c['Key for dropdown entries']]
            is_list = tokens[c['List?']]
            list_delimiter = tokens[c['Delimiter']]

            # paper: paper terms
            # mutant: mutant terms
            # mutation: mutation data
            # hidden: human curator info
            # gechange: transcriptomic data

            term_list_dict[category].append(entry_name)

            parsing_dict = {}
            parsing_dict['type'] = input_type
            parsing_dict['is_list'] = is_list == '1'
            parsing_dict['list_key'] = list_delimiter
            parsing_dict['required'] = input_required == '1'
            parsing_dict['is_dropdown'] = input_entry.upper() != 'NONE'
            parsing_dict['dropdown_key'] = input_entry
            parsing_dict['multiple'] = multiple_values == '1'
            parsing_dict['category'] = 'not used here'

            variable_types[entry_name] = input_type
            parser_rules[entry_name] = parsing_dict

    return dict(term_list_dict), variable_types, parser_rules


def parse_annotation_format_rules():

    parser_rules = {}

    annotation_terms = open(os.path.join(constants.INPUT_DIR,
                                         'settings',
                                         'Annotation Grammar_typed.txt'), 'r')

    with annotation_terms as f:

        headers = f.readline().strip().split('\t')
        c = {}

        for i, token in zip(range(0, len(headers)), headers):
            c[token] = i

        # used by MetEngDatabase to parse annotation data
        for line in f:

            tokens = line.strip().split('\t')

            mutation_name = tokens[c['Mutation']]
            mutation_type = tokens[c['Type']]
            category = tokens[c['Category']]

            # 0 = false, 1 = true
            is_list = tokens[c['List']] == '1'
            is_tuple = tokens[c['Tuple']] == '1'

            list_delimiter = tokens[c['List_Delimiter']]
            tuple_delimiter = tokens[c['Tuple_Delimiter']]
            tuple_types = tokens[c['Tuple_Types']]

            if list_delimiter.upper() == 'NONE':
                list_delimiter = ''

            if tuple_delimiter.upper() == 'NONE':
                tuple_delimiter = ''

            parsing_dict = {'type': mutation_type,
                            'is_list': is_list,
                            'is_tuple': is_tuple,
                            'list_key': list_delimiter,
                            'tuple_key': tuple_delimiter,
                            'category': category,
                            'type_formatter': tuple_types.split(',')}

            parser_rules[mutation_name] = parsing_dict

    return parser_rules


def load_tolerance_ontology():

    """
    
    Process the 'Stress Classifications.txt' file and return:
    
    tagged: dict ( str : list) where the keys are phenotypes and the value classification tags (solvent_biofuels, etc)
    ontology: tuple (str, str) where the first entry is the highest level classification, and the second entry
    the general category of the stressor (solvent_biofuels for n-butanol)
    
    :return: 
    """

    fhandle = open(os.path.join(constants.INPUT_DIR,
                                'standardization',
                                'Stress Classifications.txt'), 'r')
    lines = fhandle.readlines()
    fhandle.close()

    tags = {}
    ontology = {}

    for line in lines[1:]:

        tokens = line.strip().replace('\"', '').split('\t')

        tokens = [x.strip() for x in tokens]

        try:
            stressor = tokens[0].upper().replace('_sensitive'.upper(), '')
            tagged_info = tokens[1]
            general_category = set([x.upper() for x in tokens[2].split(',')])
            root_category = tokens[3]
        except IndexError:
            raise IndexError('Error in format of line in Stress Classifications.txt: %s' % line)

        ontology[stressor.upper()] = (root_category, general_category)
        tags[stressor.upper()] = tagged_info

        # all phenotypes given the above removal of _sensitive
        if '_sensitive'.upper() not in stressor:
            tags[(stressor + '_sensitive').upper()] = tagged_info
            ontology[(stressor + '_sensitive').upper()] = (root_category, general_category)

    return tags, ontology


def load_fasta(name):

    """
    
    Simple FASTA file format parser.
    
    :param name: 
    :return: 
    """

    gene_sequence_dict = defaultdict(list)
    current_key = ''

    with open(name, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if '>' in line:
                current_key = line
            elif current_key != '':
                # avoids in place construction of combined string object
                gene_sequence_dict[current_key].append(line)

    output = dict()
    for key in gene_sequence_dict:
        output[key] = ''.join(gene_sequence_dict[key])

    return output


def get_species_standard():

    """
    
    Returns a Standard object meant to convert species specific names from 'Species-Specific Gene Names.txt' to 
    a standardized, but not aggregated format (i.e. their specific accessions).
    
    :return: 
    """

    return Standard(names_to_standards_file=os.path.join(constants.INPUT_DIR,
                                                         'standardization',
                                                         'Species-Specific Gene Names_DB.txt'))


def get_unified_standard():

    """

    Returns a Standard object meant to convert species specific names from 'Species-Specific Gene Names.txt' and
    Sequece-Specific Gene Names.txt into MG1655 accessions.

    :return: 
    """

    species_std = Standard(names_to_standards_file=os.path.join(constants.INPUT_DIR,
                                                                'standardization',
                                                                'Species-Specific Gene Names_DB.txt'))
    unified_std = Standard(names_to_standards_file=os.path.join(constants.INPUT_DIR,
                                                                'standardization',
                                                                'Sequence-Specific Gene Names_DB.txt'))

    species_std.overwrite(unified_std)

    return species_std


def get_phenotype_standard():

    """
    
    Loads some basic compound information consisting of common name and a standardized name, returns a standard
    object for handling conversions.
    
    :return: 
    """

    return Standard(os.path.join(constants.INPUT_DIR,
                                 'standardization',
                                 'Compound Information.txt'))


def get_database(species_standard=None,
                 phenotype_standard=None,
                 get_parser=False,
                 folders=('2015', '2016', '2017', '2018')):

    """
    
    Returns a set of records consisting of Paper, Mutant, Mutation, Transcriptomic, and Annotation objects
    that represented parsed text records. It is possible to directly use these objects for analysis, but 
    you should use SQL to interface with the database instead.
    
    This function is mainly used by the database sql parsing script to generate the associated SQL database.
    
    Note that the folders flag is actually the directories read by the function to gather data for parsing,
    and it can include other directories within the database_store/ directory.
    
    :param species_standard: 
    :param phenotype_standard: 
    :param get_parser: 
    :param folders: 
    :return: 
    """

    if species_standard is None:
        species_standard = get_species_standard()

    if phenotype_standard is None:
        phenotype_standard = get_phenotype_standard()

    annotation_parser_rules = parse_annotation_format_rules()

    term_list_dict, v_types, term_parser_rules = parse_input_terms()

    parser_rules = {x: y for x,y in annotation_parser_rules.items()}
    for x in term_parser_rules:
        parser_rules[x] = term_parser_rules[x]

    data_terms = term_list_dict['paper']
    mutant_terms = term_list_dict['mutant']
    gene_terms = term_list_dict['gene']
    ge_terms = term_list_dict['gechange']

    # fields that represent the number of dependent mutants, mutations, and GEchanges
    # in each record
    num_mk = 'NumberofMutants'
    num_gk = 'NumberofMutations'
    num_ge = 'NumberofGEChanges'

    paths = []

    for folder in folders:
        # read files in database store directory
        # does not include all readable files
        # excludes Carol Gross data (not for scientific reasons, can be added if wanted)
        path = os.path.join(constants.INPUT_DIR,
                            'database_store',
                            folder)

        record_file_list = os.listdir(path)
        paths.extend([os.path.join(path, x) for x in record_file_list])

    database = ContainerClass(paths,
                              num_mk,
                              num_gk,
                              num_ge,
                              data_terms,
                              mutant_terms,
                              gene_terms,
                              ge_terms,
                              parser_rules,
                              v_types,
                              [],
                              species_standard_obj=species_standard,
                              phenotype_std=phenotype_standard)

    if get_parser:
        return database.records, parser_rules
    else:
        return database.records


def standard_mutation_formatting(mutation_type, annotation):
    """

    Generates a 'visually pleasing' version of internal Resistome mutation data depending on the mutation type
    and annotation data available.
    
    Note that visually pleasing is in the eye of the beholder.

    :param mutation_type: 
    :param annotation: 
    :return: 
    """

    if mutation_type == 'aa_snps':
        str_output = []
        for aa_array in annotation['aa_snps']:
            str_output.append(aa_array[1] + str(aa_array[0]) + aa_array[2])
        return 'AA change(s):' + ', '.join(str_output)
    elif mutation_type == 'nuc_snps':
        str_output = []
        for nuc_array in annotation['nuc_snps']:
            str_output.append(nuc_array[1] + str(nuc_array[0]) + nuc_array[2])
        return 'SNP(s):' + ', '.join(str_output)
    elif mutation_type == 'indel':
        str_output = []
        for indel_array in annotation['indel']:
            prefix = '+' if indel_array[1] > 0 else '-'
            try:
                str_output.append(
                    'indel:' + str(indel_array[0]) + '|' + prefix + str(abs(indel_array[1])) + ' bp|' +
                    indel_array[2])
            except:
                str_output.append('indel:location or size unclear')
        return ','.join(str_output)
    elif mutation_type == 'is_insertion':
        return 'IS insertion:' + annotation['is_insertion'][0]
    elif mutation_type == 'del':
        return 'deleted'
    elif mutation_type == 'intergenic':
        return 'intergenic:' + annotation['intergenic'][0] + '/' + annotation['intergenic'][1]
    elif mutation_type == 'amplified':
        return 'amplified:%iX' % annotation['amplified']
    elif mutation_type == 'duplication':

        output_str = 'duplication:'

        for duplication in annotation['duplication']:
            prefix = '+'
            output_str += str(duplication[0]) + '-' + str(duplication[1]) + '|' + prefix + str(
                duplication[2]) + ' bp|' + duplication[3]

        return output_str

    elif mutation_type == 'large_amplification':
        return 'large amp:' + annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (amplified %iX)' % \
                                                                                                  annotation[
                                                                                                      mutation_type][2]
    elif mutation_type == 'large_deletion':
        return 'large del:' + annotation[mutation_type][0] + '-' + annotation[mutation_type][1] + ' (deleted)'
    elif mutation_type == 'large_inversion':
        return 'large invert:' + annotation[mutation_type][0] + ' <=> ' + annotation[mutation_type][1] + ' (inverted)'
    elif mutation_type == 'mutated':
        return 'mutated'
    elif mutation_type == 'oe':
        return 'overexpressed'
    elif mutation_type == 'plasmid':
        return 'plasmid exp.'
    elif mutation_type == 'truncated':
        return 'truncated'
    elif mutation_type == 'rep':
        return 'repressed'
    elif mutation_type == 'integrated':
        return 'integrated'
    elif mutation_type == 'frameshift':
        return 'frameshift'
    elif mutation_type == 'rbs_tuned':
        return 'engineered RBS'
    else:
        return '(%s)' % mutation_type

