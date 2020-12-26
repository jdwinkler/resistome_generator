import os
from collections import defaultdict

current_location = os.path.realpath(__file__)
path = os.path.split(os.path.split(current_location)[0])[0]
base_dir = path

INPUT_DIR = os.path.join(base_dir, 'inputs')
OUTPUT_DIR = os.path.join(base_dir, 'output')
DATABASE_DIR = os.path.join(base_dir, 'database_store')
ERROR_DIR = os.path.join(base_dir, 'errors')

# define universal constants relating to mutation classification
# this is to standardize these definitions across the package
# this mutation terms are listed in the Terms Usage.txt file
OE_mutations = {'amplified', 'oe', 'plasmid', 'rbs_tuned', 'large_amplification'}
REP_mutations = {'rep', 'antisense'}
DEL_mutations = {'del', 'frameshift', 'large_deletion'}
RANDOM_mutations = {'aa_snps', 'nuc_snps', 'intergenic', 'mutated'}
REARRANGEMENT_mutations = {'large_inversion'}
INDEL_mutations = {'is_insertion', 'indel'}
ADD_mutations = {'sensor', 'regulated'}

SPECIES_LIST = {'w', 'rel606', 'mg1655', 'w3110', 'bw25113', 'mds42', 'bl21', 'bl21_de3_'}
MAP_SPECIES_TO_REF = {'K-12': 'mg1655', 'K12': 'mg1655'}

reference_genomes = [('mg1655', '.2', 2013, 'NC_000913.2.ptt')]

reference_genome_dict = defaultdict(dict)
select_reference_genome_filter = defaultdict(list)

for (species, assembly_name, assembly_year, reference_genome_data) in reference_genomes:

    reference_genome_dict[species][assembly_name] = dict()
    select_reference_genome_filter[species].append((assembly_year, assembly_name))

    with open(os.path.join(INPUT_DIR, 'biological_info', reference_genome_data)) as f:

        columns = f.readline()
        split_columns = columns.strip().split('\t')
        column_dict = {column: i for i, column in zip(range(0, len(split_columns)), split_columns)}

        for line in f:

            tokens = line.strip().split('\t')
            location_data = tokens[column_dict['Location']].split('..')
            (start, stop) = (location_data[0], location_data[1])
            accession = tokens[column_dict['Synonym']].upper()

            reference_genome_dict[species][assembly_name][accession] = (int(start), int(stop))

# significance threshold for statistical methods
# maximum threshold for significance, can be lowered for multi-hypothesis testing automatically
P_VALUE_THRESHOLD = 0.05

NUCLEOTIDE_BASES = {'A', 'T', 'C', 'G'}
AMINO_ACIDS = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', '*']

with open(os.path.join(INPUT_DIR, 'settings', 'Term Usage.txt'), 'r') as f:

    collation_dict = defaultdict(set)

    for line in f:
        tokens = line.strip().split('\t')
        collation_dict[tokens[0]].add(tokens[1])

    AIR = collation_dict['Air']
    CULTURE_VESSELS = collation_dict['Culture']
    EXPERIMENTAL_DESIGNS = collation_dict['DOE']
    MUTATION_EFFECTS = collation_dict['Effect']
    ENGINEERING_METHODS = collation_dict['Methods']
    SCORE_REASONS = collation_dict['Reason']
    GROWTH_PHASES = collation_dict['GrowthPhase']
    EXPRESSION_METHOD = collation_dict['GEMethod']


with open(os.path.join(INPUT_DIR, 'db_credentials', 'credentials.txt')) as f:

    DB_USERNAME = f.readline().strip().split('=')[1]
    DB_PASSWORD = f.readline().strip().split('=')[1]
    DB_NAME = f.readline().strip().split('=')[1]


def get_strain_converter(strain):
    """

    Converts strain (str) into the following:

    If Str is REL606 or B, returns REL606
    If W, return W
    Otherwise, assume the strain is K-12 and return MG1655.

    :param strain: 
    :return: 
    """

    strain = strain.upper()

    if strain in {'REL606', 'B'}:
        return 'rel606'
    if strain in {'W', 'BW25113', 'MDS42', 'MG1655', 'BL21', 'W3110'}:
        # return strain as-is since we have gene mappings
        return strain.lower()
    if strain == 'K-12' or strain == 'K12':
        return 'mg1655'
    if strain == 'BL21(DE3)':
        # ( ) are reserved chars in most DBMSs
        return 'bl21_de3_'

    # default. might be wrong in some cases.
    return 'mg1655'


def get_proper_start_stop(gene, study_year, query_species):

    if query_species not in SPECIES_LIST:
        raise AssertionError('Unknown species requested: %s' % query_species)

    if query_species not in select_reference_genome_filter:
        return None

    available_references = select_reference_genome_filter[query_species]

    reference = None

    for (year_cutoff, assembly_name) in available_references:
        if study_year <= year_cutoff:
            reference = assembly_name

    if reference is not None and gene in reference_genome_dict[query_species][reference]:
        return reference_genome_dict[query_species][reference][gene]
    else:
        return None


def convert_to_relative_position(gene_start, gene_stop, position):

    gene_length = gene_stop - gene_start

    if position is None:
        return None

    try:
        position = int(position)
    except (ValueError, TypeError):
        return None

    if position < 0 or position > gene_length:
        return position
    else:
        return position - gene_start


def convert_to_absolute_position(gene_start, gene_end, position, fix_one_indexing=True):

    """
    
    Converts all types of location or position data into zero-indexed, absolute locations. Will return None
    if the converted location would be ambiguous.
    
    :param gene_start: 
    :param gene_end: 
    :param position: 
    :param fix_one_indexing: 
    :return: 
    """

    offset_to_subtract = 1 if fix_one_indexing else 0

    if gene_end > gene_start:
        gene_length = gene_end - gene_start
    else:
        return None

    if position is None:
        return None

    try:
        position = int(position)
    except (ValueError, TypeError):
        return None

    if position <= gene_length:
        return position + gene_start - offset_to_subtract
    else:
        return position - offset_to_subtract


def complement(base):

    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'G':
        return 'C'
    if base == 'C':
        return 'G'

    raise AssertionError('Invalid base: %s for complementing' % base)
