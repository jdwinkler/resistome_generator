"""

Class implementing wrapper objects for representing Resistome text files.

All file parsing is done in ContainerClass, though the constructors of each
object class handle any particular idiosyncrasies of object; I have attempted
to minimize these as much as possible.

The object hierarchy is:

Paper
    Mutant
        Mutation
        Transcriptome

Note that there is no minimum number of objects; you can have a Paper class
with no mutants, mutants with no mutations nor transcriptomic data, etc.
However, all Mutants currently have at least one mutation or transcriptome
object associated.

"""

from resistome.utils.input_standardization import Standard


class ContainerClass:

    """
    
    Wrapper class to handle storage of Paper-Mutant-Mutation/Transcriptome objects. Duplicate records,
    as judged by DOI, are not allowed.
    
    """

    # default values for missing fields

    defaults = {'text': None,
                'float': None,
                'integer': None,
                'int': None,
                'str': None}

    def __init__(self, file_paths,
                 num_mutants_key,
                 num_mutations_key,
                 num_ge_key,
                 paper_terms,
                 mutant_terms,
                 gene_terms,
                 gene_expression_terms,
                 parsing_rules,
                 type_map,
                 known_mutations,
                 species_standard_obj=None,
                 phenotype_std=None):

        """
        
        Handles the mechanics of parsing (finds the files, passes fields to appropriate constructors.
        
        :param file_paths: Resistome text records to parse (absolute file paths)
        :param num_mutants_key: name of field representing number of mutants in a Resistome record
        :param num_mutations_key: same, but for mutations
        :param num_ge_key: same, but for gene expression data
        :param paper_terms: fields comprising a Paper object / fields in the text file
        :param mutant_terms: same, but for mutant objects
        :param gene_terms: same, but for mutation objects
        :param gene_expression_terms: same, but for gene expression terms
        :param parsing_rules: dict describing how to parse each mutation (elaborated on below)
        :param type_map: type of each field
        :param known_mutations: mutations to ignore/present in background strain
        :param species_standard_obj: Standard object to convert gene names
        :param phenotype_std: Standard object to convert phenotype names
        
        """

        # unique values of certain entries (title and doi)
        self.uniqueValues = set()

        # list of records
        self.records = []

        # no conversion if none
        if species_standard_obj is None:
            species_standard_obj = Standard()

        if phenotype_std is None:
            phenotype_std = Standard()

        duplicates = []

        for record in file_paths:

            try:
                file_data = open(record, 'r').readlines()

                # parsing proceeds hierarchic fashion
                # Paper is parsed here
                # then all mutants
                # then all mutation/transcriptome objects
                paper = Paper(record, file_data, num_mutants_key, num_mutations_key,
                              num_ge_key,
                              paper_terms,
                              mutant_terms,
                              gene_terms,
                              gene_expression_terms,
                              parsing_rules,
                              known_mutations,
                              type_map,
                              species_standard_obj, phenotype_std)

                if paper.doi not in self.uniqueValues:
                    self.records.append(paper)
                else:
                    duplicates.append((paper.doi, record))

                self.uniqueValues.add(paper.doi)
            except:
                print("Parsing failed on DOI: " + record)
                raise

        if len(duplicates) > 0:
            for (doi, filename) in duplicates:
                print('Duplicate DOI detected: %s in file: %s' % (doi, filename))
            raise AssertionError

    @staticmethod
    def parse_record(file_data):

        """
        
        Parses an entire resistome record into a dict (field name : value).
        
        :param file_data: 
        :return: 
        
        """

        backing = {}
        key_order = []

        for line in file_data:

            if "#" not in line and len(line.strip()) > 0:

                tokens = line.strip().split("=")

                key = tokens[0].strip()

                # some mislabeled fields in a few Resistome records
                if 'ResistanceUnit' in key and 'ResistanceUnits' not in key:
                    key = key.replace('ResistanceUnit', 'ResistanceUnits')

                key_order.append(key)
                value = tokens[1].strip().replace('"', '')

                # missing values are handled elsewhere
                if len(value) > 0:
                    backing[key] = value

        return backing, key_order

    @staticmethod
    def value_caster(term, value, variable_type):

        """
        
        Handles conversion of string value data into the proper type given in the annotation_parsing dict. List/tuple
        parsing is handled elsewhere.
        
        :param term: 
        :param value: 
        :param variable_type: 
        :return: 
        """

        # default values are a class constant
        if value == '':
            return ContainerClass.defaults[variable_type]

        if variable_type == 'text':
            return str(value)

        if variable_type == 'integer':
            return int(value)

        if variable_type == 'float':
            return float(value)

        raise AssertionError('Unrecognized variable type requested: %s, %s' % (term, variable_type))

    @staticmethod
    def parse_entry(key_array, backing, p_rules):

        """
        
        Converts the values in backing referenced by keys in key array using the rules described in p_rules
        (dict (key : dict : (str : parsing rule)).
        
        :param key_array: 
        :param backing: 
        :param p_rules: 
        :return: 
        """

        container = dict()

        for term, backing_key, storage_key in key_array:

            # attempts to pull value out of backing
            try:
                value = backing[backing_key].strip()
            except KeyError:
                # not included in backing/file
                value = ''
                container[term] = value

            if not p_rules[term]['is_list']:
                # parsing individual values according to type rules specified in p_rules
                container[storage_key] = ContainerClass.value_caster(term,
                                                                     value,
                                                                     p_rules[term]['type'])
            else:
                # parses each element of a list, assumes comma separator
                if value != '':
                    container[storage_key] = [
                        ContainerClass.value_caster(term, x.strip(),
                                                    p_rules[term]['type'])
                        for x in value.split(",")]
                else:
                    # empty list otherwise
                    container[storage_key] = []
        return container

    @staticmethod
    def parse_paper_entries(backing, paper_term_array, p_rules):

        """
        
        Parses all paper entries using whatever is listed in paper_term_array.
        
        :param backing: 
        :param paper_term_array: 
        :param p_rules: 
        :return: 
        """

        paper_backing = ContainerClass.parse_entry(paper_term_array, backing, p_rules)
        return paper_backing

    @staticmethod
    def parse_mutant_entries(backing, num_mutants, mutant_terms, p_rules):

        """
        
        Parses all Mutant related entries in backing. Converts mutant terms into Mutant#.term, where #
        is the number of mutant annotated in the record as deduced from num_mutants.
        
        :param backing: 
        :param num_mutants: 
        :param mutant_terms: 
        :param p_rules: 
        :return: 
        """

        key_array = []
        for i in range(1, num_mutants + 1):
            for term in mutant_terms:
                # pattern: Mutant-number.term
                key = 'Mutant' + str(i) + '.' + term
                key_array.append((term, key, (i, term)))

        mutant_backing = ContainerClass.parse_entry(key_array, backing, p_rules)

        return mutant_backing

    @staticmethod
    def parse_gene_entries(backing, mutant_number_array, mutation_count_array,
                           gene_terms, p_rules):

        """
        
        Parse gene and annotation entries (since they are not independent). Also converts gene and annotation terms
        as follows: Mutant-number_i.Mutation-number_j.term where number_i is the mutant number, and number_j is the
        mutation number.
        
        :param backing: 
        :param mutant_number_array: 
        :param mutation_count_array: 
        :param gene_terms: 
        :param p_rules: 
        :return: 
        """

        gene_mutation_positions = []
        key_array = []
        for i, j in zip(mutant_number_array, mutation_count_array):
            for term in gene_terms:

                key = 'Mutant' + str(i) + '.Mutation' + str(j) + '.' + term
                key_array.append((term, key, (i, j, term)))

                # collate for annotation parsing
                if term == 'GeneMutation':
                    gene_mutation_positions.append((i, j, term))

        container = ContainerClass.parse_entry(key_array,
                                               backing,
                                               p_rules)

        annotation_dict = ContainerClass.parse_annotations(container,
                                                           backing,
                                                           gene_mutation_positions,
                                                           p_rules)

        return container, annotation_dict

    @staticmethod
    def parse_transcript_entries(backing, mutant_number_array,
                                 expression_count_array, expression_terms,
                                 p_rules):

        """
        
        Same as parse_gene_entries, but for transcriptomic changes instead. Queries backing for 
        Mutant-number.GEChange-number.term entries.
        
        :param backing: 
        :param mutant_number_array: 
        :param expression_count_array: 
        :param expression_terms: 
        :param p_rules: 
        :return: 
        """

        key_array = []
        for i, j in zip(mutant_number_array, expression_count_array):
            for term in expression_terms:
                key = 'Mutant' + str(i) + '.GEChange' + str(j) + '.' + term
                key_array.append((term, key, (i, j, term)))

        container = ContainerClass.parse_entry(key_array, backing, p_rules)

        return container

    @staticmethod
    def parse_annotations(gene_backing, backing, tuple_keys, p_rules):

        """
        
        Parses annotations into a consistent format as described by p_rules. 
        
        :param gene_backing: 
        :param backing: 
        :param tuple_keys: 
        :param p_rules: 
        :return: 
        """

        defaults = {'str': None, 'int': float('nan')}

        def tuple_value_caster(vtuple, requested_types):

            """
            
            Helper function for parsing and casting annotation values.
            
            :param vtuple: 
            :param requested_types: 
            :return: 
            """

            items = []
            boundary = len(vtuple) if len(vtuple) > len(requested_types) else len(requested_types)
            for j in range(0, boundary):
                if j >= len(vtuple):
                    # if a value is missing, write in the default defined above
                    # generally this is pretty rare
                    items.append(defaults[requested_types[j]])
                else:
                    try:
                        if requested_types[j] == 'int':
                            # locations
                            converted_value = int(vtuple[j])
                        else:
                            # default is generally string
                            converted_value = str(vtuple[j])
                        items.append(converted_value)
                    except (ValueError, IndexError):
                        # use str version of vutple entry
                        items.append(str(vtuple[j]))
            return tuple(items)

        def apply_parsing_rules(mutation_x, value_x):

            if mutation_x not in p_rules:
                # new mutation types need to be define in InputTypes / Term Usage
                raise AssertionError('Missing rule for %s mutation annotation' % mutation_x)
            parsing_dict = p_rules[mutation_x]
            if parsing_dict['is_list'] and parsing_dict['is_tuple']:
                # these are parsed in an identical fashion
                # can have multiple tuples though (not exclusive)
                tuple_tokens = [x.strip() for x in value_x.split(parsing_dict['list_key'])]
                output = []
                for y in tuple_tokens:
                    tokens = [x.strip() for x in
                              y.split(parsing_dict['tuple_key'])]
                    tokens = [ContainerClass.value_caster(mutation_x, x,
                                                          parsing_dict['type'])
                              for x in tokens]
                    output.append(tuple_value_caster(tokens, parsing_dict['type_formatter']))

                return output

            elif parsing_dict['is_list'] or parsing_dict['is_tuple']:
                # either or but not both
                if parsing_dict['is_list']:
                    delimiter = parsing_dict['list_key']
                    tokens = [x.strip() for x in value_x.split(delimiter)]
                    tokens = [ContainerClass.value_caster(mutation_x, x,
                                                          parsing_dict['type'])
                              for x in tokens]
                    return tokens
                else:
                    delimiter = parsing_dict['tuple_key']
                    tokens = [x.strip() for x in value_x.split(delimiter)]
                    tokens = [ContainerClass.value_caster(mutation_x, x,
                                                          parsing_dict['type'])
                              for x in tokens]
                    return tuple_value_caster(tokens,
                                              parsing_dict['type_formatter'])
            else:
                return ContainerClass.value_caster(mutation_x,
                                                   value_x,
                                                   parsing_dict['type'])

        annotation_backing = dict()

        for (i, j, term) in tuple_keys:
            mutations = gene_backing[(i, j, term)]
            for mutation in mutations:
                # mutation type forms last key
                key = 'Mutant' + str(i) + '.Mutation' + str(j) + '.' + term + '.' + mutation
                try:
                    value = backing[key]
                except KeyError:
                    value = ''

                annotation_backing[
                    (i, j, term, mutation.lower())] = apply_parsing_rules(
                    mutation, value)

        return annotation_backing


class Paper:
    def __init__(self, file_name,
                 file_data,
                 num_mutants_key,
                 num_mutations_key,
                 num_ge_key,
                 paper_terms,
                 mutant_terms,
                 gene_terms,
                 expression_terms,
                 p_rules,
                 known_mutations,
                 type_map,
                 standard_obj,
                 phenotype_std):

        """
        
        Wrapper class for Paper data, and maintains arrays of mutant objects.
        
        :param file_name: 
        :param file_data: 
        :param num_mutants_key: 
        :param num_mutations_key: 
        :param num_ge_key: 
        :param paper_terms: 
        :param mutant_terms: 
        :param gene_terms: 
        :param expression_terms: 
        :param p_rules: 
        :param known_mutations: 
        :param type_map: 
        :param standard_obj: 
        :param phenotype_std: 
        """

        # backing contains the original key-value file, but with no
        # parsing of the value beyond space trimming.
        self.fileName = file_name
        self.keyOrder = []

        self.parsingRules = p_rules

        self.type_map = type_map

        backing, self.keyOrder = ContainerClass.parse_record(file_data)
        paper_backing = ContainerClass.parse_paper_entries(backing, [(x, x, x) for x in paper_terms], p_rules)
        num_mutants = int(backing.get(num_mutants_key, 0))

        mutant_backing = ContainerClass.parse_mutant_entries(backing, num_mutants, mutant_terms, p_rules)

        # reference genomes change over time
        # if a record annotates this, record it, otherwise use default latest
        self.genome_reference = backing.get('genome_reference', 'U00096.3')

        mutant_array_genes = []
        mutation_array = []
        expression_array = []

        for i in range(1, num_mutants + 1):

            if mutant_backing[(i, num_mutations_key)] is not None:
                num_mutations = int(mutant_backing[(i, num_mutations_key)])
            else:
                # if no mutational data is recorded
                num_mutations = 0

            for j in range(1, num_mutations + 1):
                mutant_array_genes.append(i)
                mutation_array.append(j)

        mutant_array_transcript = []

        for i in range(1, num_mutants + 1):

            if mutant_backing[(i, num_ge_key)] is not None:
                num_ge_changes = int(mutant_backing[(i, num_ge_key)])
            else:
                # if no mutational data is recorded
                num_ge_changes = 0

            for j in range(1, num_ge_changes + 1):
                mutant_array_transcript.append(i)
                expression_array.append(j)

        # parse backing dicts for gene/annotation data
        gene_backing, annotation_backing = ContainerClass.parse_gene_entries(
            backing,
            mutant_array_genes,
            mutation_array,
            gene_terms,
            p_rules)

        # parse backing dicts for expression data
        expression_backing = ContainerClass.parse_transcript_entries(
            backing,
            mutant_array_transcript,
            expression_array,
            expression_terms,
            p_rules)

        self.mutants_x = []

        for i in range(1, num_mutants + 1):
            mutant = Mutant(paper_backing['Title'], i,
                            mutant_backing,
                            gene_backing,
                            expression_backing,
                            annotation_backing,
                            mutant_terms,
                            gene_terms,
                            expression_terms,
                            standard_obj,
                            phenotype_std,
                            p_rules,
                            num_gk=num_mutations_key)

            self.mutants_x.append(mutant)

        for key in paper_terms:
            setattr(self, key, paper_backing[key])

    def reference_genome(self):
        # returns reference genome
        return self.genome_reference

    def get_mutant(self, i):
        # returns mutant object i (corrects for 1 indexing in records)
        return self.mutants[i - 1]

    @property
    def mutants(self):
        # list of mutant objects
        return self.mutants_x

    @property
    def doi(self):
        # paper doi
        return getattr(self, 'DOI')

    @property
    def title(self):
        # paper title
        return getattr(self, 'Title')

    @property
    def year(self):
        # year
        return getattr(self, 'Year')

    @property
    def rgroup(self):
        # name of PI
        return getattr(self, 'ResearchGroup')

    @property
    def journal(self):
        # journal published in
        return getattr(self, 'Journal')

    @property
    def difficulty(self):
        # project difficulty score (1: not challenging, 5: extremely difficult)
        return getattr(self, 'project_score')

    @property
    def difficulty_reason(self):
        # reason for project difficulty score
        return getattr(self, 'score_reason')

    @property
    def design_method(self):
        # design method
        return getattr(self, 'project_doe')

    @property
    def categories(self):
        # tags indicating category (metal_ions, solvents_biofuels, etc)
        return getattr(self, 'Tags')

    @property
    def total_designs(self):
        # total number of designs
        return getattr(self, 'total')

    @property
    def comments(self):
        # free form text comments
        return getattr(self, 'pComments')


class Mutant:
    def __init__(self, title,
                 mutant_number,
                 mutant_backing,
                 gene_backing,
                 expression_backing,
                 annotation_backing,
                 mutant_fields,
                 gene_fields,
                 expression_fields,
                 standard_obj,
                 phenotype_std,
                 p_rules,
                 num_gk='NumberofMutations',
                 num_ge='NumberofGEChanges',
                 ignore_values=False):

        """
        
        Wrapper class for Mutant object data, and maintains records for mutation, annotation, and gene expression
        objects.
        
        :param title: 
        :param mutant_number: 
        :param mutant_backing: 
        :param gene_backing: 
        :param expression_backing: 
        :param annotation_backing: 
        :param mutant_fields: 
        :param gene_fields: 
        :param expression_fields: 
        :param standard_obj: 
        :param phenotype_std: 
        :param p_rules: 
        :param num_gk: 
        :param num_ge: 
        :param ignore_values: 
        """

        self.id = mutant_number
        self.backing = {}
        self.mutations = []
        self.expression = []
        all_affected_genes = set()
        ancillary_genes = set()

        # copies mutant i information into it's own data structure
        for term in mutant_fields:

            if ignore_values and (self.id, term) not in mutant_backing:
                # default empty value
                self.backing[term] = ''
            else:
                if term == 'TolerancePhenotype' or term == 'AntagonisticPhenotype':
                    # sensitive phenotypes should be in AntagonisticPhenotype, and are added there later
                    # remove just removes the text '_sensitive' when checking the Standard object for this conversion
                    # it is added back on to the converted name if present
                    try:
                        self.backing[term] = [phenotype_std.convert(x, 'None', remove='_sensitive').upper()
                                              for x in mutant_backing[(self.id, term)]]
                    except TypeError:
                        raise TypeError('Error when trying to convert Tolerance/AntagonisticPhenotype field in %s'
                                        % title)
                else:
                    self.backing[term] = mutant_backing[(self.id, term)]

        try:
            num_mutations = int(self.backing[num_gk])
        except (ValueError, TypeError):
            num_mutations = 0

        for j in range(1, num_mutations + 1):

            # filter out placeholder garbage here before making mutation object
            if ignore_values and ('GeneName' or 'GeneSource' not in self.backing):
                continue

            # you cannot have a gene named nonE from an organism called NONE
            if not (gene_backing[(self.id, j, 'GeneName')].upper() == 'NONE' or
                    gene_backing[(self.id, j, 'GeneSource')].upper() == 'NONE'):
                mutation_obj = Mutation(self.id,
                                        j,
                                        self.backing['Species'],
                                        self.backing['Subspecies'],
                                        gene_backing,
                                        annotation_backing,
                                        gene_fields,
                                        standard_obj,
                                        p_rules,
                                        ignore_missing_values=ignore_values)
                self.mutations.append(mutation_obj)

                # list of all genes affected in this mutant
                all_affected_genes.update(mutation_obj.all_affected_genes(p_rules))

                # gets other genes affected by the mutation (generally large deletions/inversion/amplifications)
                ancillary_genes.update(mutation_obj.ancillary_genes(p_rules))
        try:
            num_ge_changes = int(self.backing[num_ge])
        except (ValueError, TypeError):
            num_ge_changes = 0

        for j in range(1, num_ge_changes + 1):

            if ignore_values and ('GeneName' or 'GeneSource' not in self.backing):
                continue

            if not (expression_backing[(self.id, j, 'GeneName')].upper() == 'NONE' or
                    expression_backing[(self.id, j, 'GeneSource')].upper() == 'NONE'):
                expression_obj = Transcriptome(self.id,
                                               j,
                                               self.backing['Species'],
                                               self.backing['Subspecies'],
                                               expression_backing,
                                               expression_fields,
                                               standard_obj)

                self.expression.append(expression_obj)

        # define attributes for object to carry around
        for key in self.backing:
            if key == 'ResistanceUnit':
                # fixed in most records?
                setattr(self, '_' + key + 's', self.backing[key])
            else:
                setattr(self, '_' + key, self.backing[key])

        if self.name is None:
            # no name
            self._Name = 'mutant' + str(self.id)

        # fixes sensitive phenotypes in tolerances
        # should be in antagonistic phenotypes
        tolerances = getattr(self, '_TolerancePhenotype')

        antagonistic = self.backing['AntagonisticPhenotype']

        # remove any phenoypes from tolerances if they have _sensitive tags
        # move to antagonistic
        # store for retrival in new fields
        antagonistic.extend(
            [x.replace('_sensitive'.upper(), '') for x in tolerances if
             x.find('_sensitive'.upper()) > -1])

        setattr(self, '_sensitive_phenotypes', list(set(antagonistic)))
        setattr(self, '_resistant_phenotypes',
                [x for x in tolerances if x.find('_sensitive'.upper()) == -1])

        # reset number of mutations to match the length of mutations array
        setattr(self, num_gk, len(self.mutations))

        setattr(self, '_affected_genes', all_affected_genes)

        setattr(self, '_ancillary_genes', ancillary_genes)

        self.id = self.id - 1

    @property
    def affected_genes(self):
        # genes listed as being mutated
        return getattr(self, '_affected_genes')

    @property
    def ancillary_genes(self):
        # genes listed as being affected by the mutations but not explicitly labeled as such
        return getattr(self, '_ancillary_genes')

    def mutations(self):
        # list of mutation objects
        return self.mutations

    def expression_changes(self):
        # list of expression objects
        return self.expression

    @property
    def species(self):
        return getattr(self, '_Species')

    @property
    def strain(self):
        return getattr(self, '_Subspecies')

    @property
    def culture_vessel(self):
        return getattr(self, '_CultureSystem')

    @property
    def oxygenation(self):
        # aerobic
        # microaerobic
        # anaerobic
        return getattr(self, '_Oxygen')

    @property
    def medium(self):
        # (medium descriptor, [any supplements list])
        return getattr(self, '_Medium'), getattr(self, '_Supplements')

    @property
    def culture_volume(self):
        try:
            return float(getattr(self, '_cvolume'))
        except (ValueError, TypeError):
            return getattr(self, '_cvolume')

    @property
    def liquid_volume(self):

        try:
            return float(getattr(self, '_fvolume'))
        except (ValueError, TypeError):
            return getattr(self, '_fvolume')

    @property
    def substrates(self):
        return getattr(self, '_CarbonSource')

    @property
    def pH(self):
        try:
            return float(getattr(self, '_pH'))
        except (ValueError, TypeError):
            return getattr(self, '_pH')

    @property
    def temperature(self):
        try:
            return float(getattr(self, '_Temperature'))
        except (ValueError, TypeError):
            return getattr(self, '_Temperature')

    @property
    def rpm(self):
        try:
            return float(getattr(self, '_Rotation'))
        except (ValueError, TypeError):
            return getattr(self, '_Rotation')

    @property
    def name(self):
        return getattr(self, '_Name')

    @property
    def methods(self):
        # method used to generate mutant
        return getattr(self, '_Method')

    @property
    def sensitive_phenotypes(self):
        # list of sensitive phenotypes
        return getattr(self, '_sensitive_phenotypes')

    @property
    def resistant_phenotypes(self):
        # list of resistance phenotypes
        return getattr(self, '_resistant_phenotypes')

    @property
    def resistance_level(self):
        # some numeric stress level
        try:
            return float(getattr(self, '_ResistanceLevel'))
        except (ValueError, TypeError):
            return getattr(self, '_ResistanceLevel')

    @property
    def resistance_units(self):
        # units of resistance
        return getattr(self, '_ResistanceUnits')

    @property
    def metabolite(self):
        # target molecule (if you are producing something)
        # this is included mainly to allow reading LASER records, which expect this field
        return getattr(self, '_TargetMolecule')

    @property
    def fold_improvement(self):
        try:
            return float(getattr(self, '_FoldImprovement'))
        except (ValueError, TypeError):
            return getattr(self, '_FoldImprovement')

    @property
    def initial_fitness(self):
        try:
            return float(getattr(self, '_initial_fitness'))
        except (ValueError, TypeError):
            return getattr(self, '_initial_fitness')

    @property
    def final_fitness(self):
        try:
            return float(getattr(self, '_final_fitness'))
        except (ValueError, TypeError):
            return getattr(self, '_final_fitness')

    @property
    def fitness_unit(self):
        return getattr(self, '_Fitness unit')

    @property
    def comments(self):
        # free form text comment
        return getattr(self, '_mComments')


class Mutation:
    def __init__(self, i, j, host, strain, gene_backing, a_dict,
                 gene_fields,
                 standard_obj, p_rules, ignore_missing_values=False):

        """
        
        Wrapper class for mutation objects, plus annotation information.
        
        :param i: 
        :param j: 
        :param host: 
        :param strain: 
        :param gene_backing: 
        :param a_dict: annotation dict
        :param gene_fields: 
        :param standard_obj: 
        :param p_rules: 
        :param ignore_missing_values: 
        """

        self.backing = {}
        self.id = j

        if standard_obj.use_species:
            species = strain
        else:
            species = None

        self.annotation_backing = {}

        for term in gene_fields:

            if ignore_missing_values and (i, j, term) not in gene_backing:
                self.backing[term] = ''
            else:
                self.backing[term] = gene_backing[(i, j, term)]

        for key in self.backing:

            if key == 'GeneName':
                setattr(self, key, standard_obj.convert(self.backing[key].replace('\'', ''),
                                                        self.backing['GeneSource'],
                                                        species))
            else:
                setattr(self, key, self.backing[key])

        for mutation in self.changes:
            if p_rules[mutation]['category'] == 'gene':
                # mutation category
                if isinstance(a_dict[(i, j, 'GeneMutation', mutation)], list) \
                        or isinstance(a_dict[(i, j, 'GeneMutation', mutation)], tuple):
                    temp_array = []
                    for x in a_dict[(i, j, 'GeneMutation', mutation)]:
                        if isinstance(x, str):
                            # convert all genes listed in mutation
                            # for larger multi-gene mutations
                            temp_array.append(standard_obj.convert(x.replace('\'', ''), host, species))
                        else:
                            # ignore for now
                            temp_array.append(x)
                    self.annotation_backing[mutation] = temp_array
                else:
                    # convert mutation name
                    self.annotation_backing[mutation] = standard_obj.convert(
                        a_dict[(i, j, 'GeneMutation', mutation)], host,
                        species)
            else:
                self.annotation_backing[mutation] = a_dict[(i, j, 'GeneMutation', mutation)]

        setattr(self, 'modified_genes', self.all_affected_genes(p_rules))

    def ancillary_genes(self, p_rules):

        """
        
        Returns all genes that are affected by this mutation, even if they are not labeled as such.
        
        For example, if rph is labeled as having an intergenic mutation between rph-pyrE, this will 
        return [pyrE].
        
        :param p_rules: 
        :return: 
        """

        all_affected_genes = self.all_affected_genes(p_rules)

        all_affected_genes.remove(self.name)

        return all_affected_genes

    def all_affected_genes(self, p_rules):

        """
        
        Returns a list of all affected genes for this mutation object.
        
        :param p_rules: 
        :return: 
        """

        output = set()
        output.add(self.name)

        for mutation in self.changes:

            if p_rules[mutation]['category'] == 'gene':

                if p_rules[mutation]['is_list'] or p_rules[mutation]['is_tuple']:
                    total_len = len(self.annotation(mutation))
                    counter = 0
                    for x in self.annotation(mutation):
                        # last entry of large_amplification is the amount of amplification
                        # not a gene name
                        if mutation == 'large_amplification' and counter == total_len - 1:
                            break
                        output.add(x)
                        counter += 1
                else:
                    output.add(self.annotation(mutation))

        return output

    def contains(self, mutation):

        return mutation in self.annotation_backing

    def annotation(self, mutation):

        if mutation in self.annotation_backing:
            return self.annotation_backing[mutation]
        else:
            raise AssertionError('Mutation not present in this mutation object: %s' % mutation)

    @property
    def source(self):
        # source strain for the gene, string
        return getattr(self, 'GeneSource')

    @property
    def is_original(self):
        # returns a boolean describing whether or not this gene has been
        # mutated de novo in this design or is common knowledge
        result = getattr(self, 'GeneOriginal')
        return result == 'yes'

    @property
    def effects(self):
        # strings describing the effect of all gene mutations together
        return getattr(self, 'GeneEffect')

    @property
    def changes(self):
        # list of mutations affecting the gene, see Term Usage.txt for
        # allowed entries (not enforced at the coding level)
        # but please don't make new mutation types without a good reason...
        return getattr(self, 'GeneMutation')

    @property
    def nicknames(self):
        # alternative names for the gene in question
        return getattr(self, 'GeneNickname')

    @property
    def confidence(self):
        # confidence that the mutation has an effect on the phenotype(s) of interest
        # very hard to define
        return getattr(self, 'EffectConfidence')

    @property
    def name(self):
        # gene name
        return getattr(self, 'GeneName')


class Transcriptome:
    def __init__(self, i, j, host, strain, ge_backing, ge_fields, standard_obj,
                 ignore_missing_values=False):

        """
        
        Similar to Mutation wrapper object, but does not incorporate any subannotation data.
        
        :param i: 
        :param j: 
        :param host: 
        :param strain: 
        :param ge_backing: 
        :param ge_fields: 
        :param standard_obj: 
        :param ignore_missing_values: 
        """

        self.backing = {}
        self.id = j

        for term in ge_fields:

            if ignore_missing_values and (i, j, term) not in ge_backing:
                self.backing[term] = ''
            else:
                self.backing[term] = ge_backing[(i, j, term)]

        for key in self.backing:

            if key == 'GeneName':
                setattr(self, key, standard_obj.convert(self.backing[key], host,
                                                        strain=strain))
            else:
                setattr(self, key, self.backing[key])

    @property
    def name(self):
        # gene name
        return getattr(self, 'GeneName')

    @property
    def gechange(self):
        # overexpressed or underexpressed
        return getattr(self, 'GEChange')

    @property
    def pvalue(self):
        # p-value (adjusted)
        return getattr(self, 'Pvalue')

    @property
    def stattest(self):
        # statistical test used
        return getattr(self, 'StatisticalTest')

    @property
    def foldchange(self):
        # fold change value
        return getattr(self, 'FoldChange')

    @property
    def exposure(self):
        # exposure time when RNA was taken (hrs)
        return getattr(self, 'ExposureTime')

    @property
    def growthphase(self):
        # growth phase when RNA was taken
        return getattr(self, 'GrowthPhase')

    @property
    def gemethod(self):
        # method = RNAseq, microarray, qPCR, etc
        return getattr(self, 'GEMethod')

    @property
    def accession(self):
        # accession number and repository, ie GEO GSE32131
        return getattr(self, 'AccessionNo')

    @property
    def genesource(self):
        # gene source, like Escherichia coli, or Acinetobacter baumannii, etc
        return getattr(self, 'GeneSource')

    @property
    def stress_amount(self):
        # degree of stress
        # usually numeric?
        return getattr(self, 'StressAmount')

    @property
    def stress_units(self):
        # units of stress
        return getattr(self, 'StressUnits')
