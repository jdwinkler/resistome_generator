import sys
from collections import defaultdict

# terms used for the paper object

paper_terms = ['Title',
               'DOI',
               'Year',
               'ResearchGroup',
               'Journal',
               'project_score',
               'score_reason',
               'project_doe',
               'Tags',
               'pComments',
               'total',
               'NumberofMutants',
               'GEQC']

# terms used for gene objects

gene_terms = ['GeneName',
              'GeneNickname',
              'GeneMutation',
              'GeneEffect',
              'EffectConfidence',
              'GeneSource',
              'GeneOriginal',
              'EcNumber',
              'gComments']

# terms used for gene expression objects

ge_terms = ['GeneName',
            'GEChange',
            'Pvalue',
            'FoldChange',
            'ExposureTime',
            'GrowthPhase',
            'GEMethod',
            'AccessionNo',
            'GeneSource',
            'StressAmount',
            'StressUnits',
            'StatisticalTest']

# terms used to describe mutants

mutant_terms = ['Species',
                'Subspecies',
                'CultureSystem',
                'Oxygen',
                'Medium',
                'CarbonSource',
                'Supplements',
                'cvolume',
                'fvolume',
                'pH',
                'Temperature',
                'Rotation',
                'Name',
                'Method',
                'TolerancePhenotype',
                'AntagonisticPhenotype',
                'ResistanceLevel',
                'ResistanceUnit',
                'TargetMolecule',
                'Objective',
                'Pathway',
                'FoldImprovement',
                'initial_fitness',
                'final_fitness',
                'Fitness unit',
                'mComments',
                'NumberofMutations',
                'NumberofGEChanges']

NUCLEOTIDE_BASES = {'A', 'T', 'C', 'G'}
AMINO_ACIDS = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'S', 'T', 'C', 'M', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', '*']


class MutationHandler():

    def process_genomic_study(self, filename):

        """
        
        Handles the mechanics of loading the genetic data, parsing it, and then constructing the requested
        output dicts.
        
        Override by subclassing.
        
        :param filename: 
        :return: 
        """

        fhandle = open(filename, 'r')
        lines = fhandle.readlines()
        fhandle.close()

        delimiter = '\t'

        design_dict = defaultdict(list)

        column_headers = lines[0].strip().replace('\x00', '').replace('\xff', '').replace('\xfe', '').split(delimiter)

        c = {}
        pos = {}
        for col, i in zip(column_headers, range(0, len(column_headers))):
            c[col] = i
            pos[i] = col

        output_dict = defaultdict(list)

        for line in lines[1:]:
            line = line.strip().replace('\x00', '')

            tokens = line.split(delimiter)

            design_key = tokens[0]

            gene = tokens[c['Gene']]
            mutation_type = tokens[c['Mutation']]
            additional_info = tokens[c['Annotation']]
            location = ''

            (gene_x, laser_types, annotations) = self.convert_to_laser_mutations(gene,
                                                                            location,
                                                                            mutation_type,
                                                                            additional_info)

            design_dict[design_key].append((gene_x, laser_types, annotations))

        output_design_dict = defaultdict(list)

        for design in design_dict:

            compressed_types = defaultdict(set)
            compressed_annotations = defaultdict(dict)
            mutation_list = design_dict[design]
            for (gene, laser_types, annotations) in mutation_list:
                for type in laser_types:
                    compressed_types[gene].add(type)
                for laser_type in annotations:
                    if (laser_type not in compressed_annotations[gene]):
                        compressed_annotations[gene][laser_type] = set()

                    compressed_annotations[gene][laser_type].add(annotations[laser_type])

            for gene in compressed_types:
                output_design_dict[design].append((gene, compressed_types[gene], compressed_annotations[gene]))

        return output_design_dict, output_dict

    def process_expression_study(self, filename, threshold=0.05):

        """
        
        Parses any expression data in the input files. Threshold refers to the p-value threshold for data inclusion.
        Override by subclassing.
        
        :param filename: 
        :param threshold: 
        :return: 
        """

        fhandle = open(filename, 'r')
        lines = fhandle.readlines()
        fhandle.close()

        delimiter = '\t'

        design_dict = defaultdict(list)

        column_headers = lines[0].strip().split(delimiter)

        c = {}
        pos = {}
        for col, i in zip(column_headers, range(0, len(column_headers))):
            c[col] = i
            pos[i] = col

        for line in lines[1:]:

            line = line.strip()
            tokens = line.split(delimiter)
            gene_id = tokens[c['ID']]

            for i in range(c['CV101'], c['CV116'] + 1):

                try:

                    design_key = pos[i]
                    log2_value = float(tokens[i])
                    pvalue = ''
                    design_dict[design_key].append((gene_id, str(log2_value), pvalue))

                except:
                    continue

        return design_dict

    def convert_to_laser_mutations(self, gene, location, mutation_type, annotations):

        laser_types = self.mutation_type_converter(gene,
                                                      mutation_type,
                                                      annotations)

        (gene, laser_annotations) = self.generate_laser_annotation(laser_types,
                                                                      gene,
                                                                      location,
                                                                      mutation_type,
                                                                      annotations)

        return gene, laser_types, laser_annotations

    def mutation_type_converter(self, gene, mutation, annotations):
        """
    
        Maps a mutation given in a paper to an internal laser type. This function is more for capturing any implicit
        mutational types that are not mentioned in an input dataset.
    
        For example, intergenic base changes are often indicates as 'ackA/ptA A2927097T' without being labeled as intergenic
        so for gene names containing slashes / , an 'integenic' effect is added to the explicitly labeled nuc_snps.
    
        These effects will vary for every paper, so you will need to change what is below in your subclass.
    
        :param gene: gene name (unprocessed)
        :param mutation: mutation type given in input dataset
        :param annotations: format depends on paper/dataset
        :return: 
        """

        effects = []

        effects.extend(mutation.split(','))

        # CHANGE BELOW AS NEEDED
        if '/' in gene:
            effects.append('intergenic')

        if 'IS' in annotations[2]:
            effects.append('is_insertion')

        if 'nuc_snps' in mutation:
            effects.append('aa_snps')

        if len(effects) == 0:
            raise KeyError('No valid LASER mapping for mutation type: %s' % gene)
        else:
            return effects

    def generate_laser_annotation(self, laser_types, gene, location, mutation_type, provided_annotations):

        """
        
        Does the heavy lifting of converting a mutation type, annotation, gene, and location into an actual laser
        annotation.
        
        Override by subclassing to add customizations or replace
        
        :param laser_types: output of mutation_type_converter
        :param gene: gene name
        :param location: absolute location if available
        :param mutation_type: type of mutation
        :param provided_annotations: annotation data
        :return: 
        """

        annotations = {}

        gene_to_return = gene

        for laser_type in laser_types:

            if laser_type == 'aa_snps':

                # format: WT residue - position - mutated_residue

                residue = provided_annotations[2]
                residue = residue[residue.find('(') + 1:residue.find(')')]
                aa_change = provided_annotations[3].split(',')
                annotations[laser_type] = check_base_or_residue(aa_converter(aa_change[0]), 'aa') \
                                          + residue + \
                                          check_base_or_residue(aa_converter(aa_change[-1]), 'aa')

            elif laser_type == 'nuc_snps':

                # format: WT base - position - Mutated base
                # position is absolute if possible

                original_nt = check_base_or_residue(provided_annotations[0], 'nucleotide')
                changed_nt = check_base_or_residue(provided_annotations[1], 'nucleotide')

                annotations[laser_type] = original_nt + location + changed_nt

            elif laser_type == 'del':

                # format: nothing

                annotations[laser_type] = ''

            elif laser_type == 'oe':

                # format: promoter name used to drive overexpression

                annotations['oe'] = 'trc'

            elif laser_type == 'plasmid':

                # format: 'high', 'medium', 'low' for copy number

                annotations['plasmid'] = 'High'

            elif laser_type == 'indel':

                # format: position | size of indel | position type = {'absolute', 'relative' to gene start}

                size = provided_annotations[2]

                prefix = '+'

                if 'IS' in size:
                    loc_tok = location.split('-')
                    size = int(loc_tok[1]) - int(loc_tok[0])
                    size = str(size)

                if ',' in provided_annotations[3]:
                    prefix = '-'

                if '-' in location:
                    indel_location = location.split('-')[0]
                else:
                    indel_location = location

                annotations[laser_type] = indel_location + '|' + prefix + size + '|absolute'

            elif laser_type == 'intergenic':

                # format : gene1 | gene2

                genes = gene.split('/')
                gene1 = genes[0]
                gene2 = genes[1]
                gene_to_return = gene1
                annotations[laser_type] = gene1 + '|' + gene2

            elif laser_type == 'is_insertion':

                # format : IS name

                is_type = provided_annotations[3]
                annotations[laser_type] = is_type

            elif laser_type == 'large_deletion':

                # format: first gene | last gene

                annotations[laser_type] = provided_annotations

            elif laser_type == 'large_inversion':

                # format: first gene | last gene

                genes = gene.split('-')
                gene1 = genes[0]
                gene2 = genes[1]
                gene_to_return = gene2
                annotations[laser_type] = gene2 + '|' + gene1

            elif laser_type == 'large_amplification':

                # format: first gene | last gene | copy number

                genes = gene.split('-')
                gene1 = genes[0]
                gene2 = genes[1]
                gene_to_return = gene2
                annotations[laser_type] = gene1 + '|' + gene2 + '|' + provided_annotations

            else:
                raise KeyError('Invalid LASER type %s passed to generate_laser_annotation_function' % laser_type)

        return (gene_to_return, annotations)


def check_base_or_residue(character, sequence_type):

    character = character.upper()

    if sequence_type not in {'aa', 'nucleotide'}:
        raise AssertionError('Invalid sequence type: %s' % sequence_type)

    if sequence_type == 'aa' and character not in AMINO_ACIDS:
        raise ValueError('Invalid AA: %s' % character)

    if sequence_type == 'nucleotide' and character not in NUCLEOTIDE_BASES:
        raise ValueError('Invalid base: %s' % character)

    return character


def aa_converter(input_aa):

    """
    Ala     A       Alanine
    Arg     R       Arginine
    Asn     N       Asparagine
    Asp     D       Aspartic acid (Aspartate)
    Cys     C       Cysteine
    Gln     Q       Glutamine
    Glu     E       Glutamic acid (Glutamate)
    Gly     G       Glycine
    His     H       Histidine
    Ile     I       Isoleucine
    Leu     L       Leucine
    Lys     K       Lysine
    Met     M       Methionine
    Phe     F       Phenylalanine
    Pro     P       Proline
    Ser     S       Serine
    Thr     T       Threonine
    Trp     W       Tryptophan
    Tyr     Y       Tyrosine
    Val     V       Valine
    Asx     B       Aspartic acid or Asparagine
    Glx     Z       Glutamine or Glutamic acid.
    Xaa     X       Any amino acid.
    TERM            termination codon
    """

    aa_code = {}
    aa_code['ALA'] = 'A'
    aa_code['ARG'] = 'R'
    aa_code['ASN'] = 'N'
    aa_code['CYS'] = 'C'
    aa_code['GLN'] = 'Q'
    aa_code['GLY'] = 'G'
    aa_code['GLU'] = 'E'
    aa_code['HIS'] = 'H'
    aa_code['ILE'] = 'I'
    aa_code['LEU'] = 'L'
    aa_code['LYS'] = 'K'
    aa_code['MET'] = 'M'
    aa_code['PHE'] = 'F'
    aa_code['PRO'] = 'P'
    aa_code['SER'] = 'S'
    aa_code['THR'] = 'T'
    aa_code['TRP'] = 'W'
    aa_code['TYR'] = 'Y'
    aa_code['VAL'] = 'V'
    aa_code['ASP'] = 'D'
    aa_code['GLX'] = 'Z'
    aa_code['XAA'] = 'X'
    aa_code['STOP'] = '*'
    aa_code['*'] = '*'

    try:
        return aa_code[input_aa.upper()]
    except KeyError:
        return input_aa.upper()


def record_to_text(data_dict, mutant_no, gene_no, terms, default_value, use_mutant=False, use_gene=False,
                   prefix='Mutation'):

    """
    
    Converts data into the Resistome format.
    
    :param data_dict: dict of (term : term value) data 
    :param mutant_no: mutant number to use
    :param gene_no: mutation number to use (also doubles as expression change number)
    :param terms: terms to pull out of data_dct
    :param default_value: default value if a term is missing from data_dict
    :param use_mutant: use mutant prefix for keys
    :param use_gene: use gene or gene expression prefix
    :param prefix: prefix to use if use_gene is true.
    :return: 
    """

    result = []

    key = ''

    if use_mutant:
        key = key + 'Mutant' + str(mutant_no) + '.'

    if use_gene:
        key = key + prefix + str(gene_no) + '.'

    for term in terms:
        final_key = key + term

        if term in data_dict:

            if isinstance(data_dict[term], list) or isinstance(data_dict[term], set):
                value = ','.join(data_dict[term])
            else:
                value = data_dict[term]

            # format is key = "<value>"

            result.append(final_key + ' = \"' + value + '\"')
        else:
            result.append(final_key + ' = \"' + default_value + '\"')

    return result

handler = MutationHandler()

designs, selection_dict = handler.process_genomic_study('genetic_data.txt')
expression_data = handler.process_expression_study('expression_data.txt')

# defined mutations reconstructed for manual verification
output_file = []

mutant_counter = 1

design_keys = list(designs.keys())
design_keys.extend(expression_data.keys())
design_keys = set(design_keys)

paper_info_dict = {}
mutant_info_dict = {}

paper_info_dict[
    'Title'] = 'E Unibus Plurum: Genomic Analysis of an Experimentally Evolved Polymorphism in Escherichia coli'
paper_info_dict['DOI'] = '10.1371/journal.pgen.1000713'
paper_info_dict['Year'] = '2009'
paper_info_dict['ResearchGroup'] = 'Frank Rosenzweig'
paper_info_dict['Journal'] = 'PLoS Genetics'
paper_info_dict['project_score'] = '2'
paper_info_dict['score_reason'] = ['scale']
paper_info_dict['project_doe'] = ['random']
paper_info_dict['Tags'] = ['general_growth']
paper_info_dict['NumberofMutants'] = str(len(design_keys))
paper_info_dict['GEQC'] = '2'

output_file.extend(record_to_text(paper_info_dict, 0, 0, paper_terms, ''))

mutant_info_dict['Species'] = 'Escherichia coli'
mutant_info_dict['Subspecies'] = 'JA122'
mutant_info_dict['Medium'] = 'Davis'
mutant_info_dict['Supplements'] = ['']
mutant_info_dict['Temperature'] = '30'
mutant_info_dict['Method'] = ['evolution']
mutant_info_dict['CarbonSource'] = 'Glucose'
mutant_info_dict['Pathway'] = ''
mutant_info_dict['mComments'] = ''
mutant_info_dict['Rotation'] = ''
mutant_info_dict['CultureSystem'] = 'chemostat'
mutant_info_dict['cvolume'] = '120'
mutant_info_dict['fvolume'] = ''
mutant_info_dict['Oxygen'] = 'aerobic'
mutant_info_dict['ResistanceLevel'] = '0.0125'
mutant_info_dict['ResistanceUnit'] = '%w/v'
mutant_info_dict['pH'] = '7'
mutant_info_dict['TolerancePhenotype'] = ['Glucose LIMITATION']
mutant_info_dict['AntagonisticPhenotype'] = ['']

output_file.append('SubmitterName = \"James Winkler\"')
output_file.append('SubmitterEmail = \"james.winkler@gmail.com\"')

for design_key in design_keys:

    design = designs.get(design_key, [])
    number_of_mutations = len(design)
    mutant_info_dict['Name'] = design_key
    mutant_info_dict['NumberofMutations'] = str(number_of_mutations)
    mutant_info_dict['NumberofGEChanges'] = str(len(expression_data.get(design_key, [])))

    # can added detailed phenotype/fitness/growth rate data
    # you need to overwrite the process genomic study method and then extract that from output dict manually

    output_file.extend(record_to_text(mutant_info_dict, mutant_counter, 0, mutant_terms, '', use_mutant=True))

    gene_counter = 1

    for (gene, laser_mutations, annotations) in design:

        gene_info_dict = {}

        # these are just default entries, change if desired

        gene_info_dict['GeneName'] = gene
        gene_info_dict['GeneEffect'] = 'metabolic_disruption'
        gene_info_dict['GeneSource'] = 'Escherichia coli'
        gene_info_dict['GeneOriginal'] = 'yes'
        gene_info_dict['GeneMutation'] = laser_mutations
        gene_info_dict['EffectConfidence'] = '3'

        output_file.extend(record_to_text(gene_info_dict, mutant_counter, gene_counter, gene_terms, '', use_mutant=True,
                                          use_gene=True))

        annotation_info_dict = {}
        terms = []

        for mut in laser_mutations:
            annotation_info_dict['GeneMutation.' + mut] = annotations[mut]
            terms.append('GeneMutation.' + mut)

        output_file.extend(
            record_to_text(annotation_info_dict, mutant_counter, gene_counter, terms, '', use_mutant=True,
                           use_gene=True))

        gene_counter += 1

    ge_counter = 1

    for (gene, fold_change, p_value) in expression_data.get(design_key, []):
        ge_dict = {'GeneName': gene,
                   'FoldChange': fold_change,
                   'GEChange': 'underexpressed' if float(fold_change) < 0 else 'overexpressed',
                   'Pvalue': p_value,
                   'ExposureTime': '',
                   'GrowthPhase': 'exponential',
                   'GEMethod': 'Microarray',
                   'AccessionNo': '',
                   'GeneSource': 'Escherichia coli JA122',
                   'StressAmount': '',
                   'StressUnits': '',
                   'StatisticalTest': ''}

        output_file.extend(record_to_text(ge_dict,
                                          mutant_counter,
                                          ge_counter,
                                          ge_terms,
                                          '',
                                          use_mutant=True,
                                          use_gene=True,
                                          prefix='GEChange'))

        ge_counter += 1

    mutant_counter += 1

fhandle = open('Record_Rozenweig_Glucose_Limitation.txt', 'w')

for line in output_file:
    fhandle.write(line + '\n')
fhandle.close()
