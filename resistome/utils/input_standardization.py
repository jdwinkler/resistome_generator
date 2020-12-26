from collections import defaultdict
from resistome import constants

__author__ = 'jdwinkler'


class Standard:

    """
    
    Automates the species specific/context specific conversion of identifiers to internal Resistome standards.
    
    """

    def __init__(self, names_to_standards_file=None):

        """
        
        Constructs Standard object. Set names_to_standards = None for an empty object (i.e. you intend to use it as 
        a dummy converter).
        
        :param names_to_standards_file: 
        """
        self.use_species = True

        if names_to_standards_file is not None:
            self.standard, self.reversion = self.parse_standard(names_to_standards_file)
        else:
            self.standard = {}

    def parse_standard(self, filename):

        """
        
        Parses a tab delimited standardization file. Format:
        
        Identifier
        Species
        Strain
        Standrd identifier 
        
        :param filename: str, path to standardization file
        :return: None
        """

        try:

            fhandle = open(filename, 'r')
            lines = fhandle.readlines()
            fhandle.close()

            reverse = defaultdict(list)

            name_dict = {}

            seen_before = set()
            exceptions = set()
            identical_checker = dict()

            for line in lines[1:]:

                tokens = line.strip().split('\t')

                tokens = [x.strip().upper() for x in tokens]

                gene_colloquial_name = tokens[0]
                species = tokens[1]
                strain = tokens[2]
                mg1655_accession = tokens[3]

                if (gene_colloquial_name, strain) in seen_before and identical_checker[(gene_colloquial_name, strain)] != mg1655_accession:
                    exceptions.add((gene_colloquial_name, strain))

                seen_before.add((gene_colloquial_name, strain))
                identical_checker[(gene_colloquial_name, strain)] = mg1655_accession

                # these columns are defined in Species Specific Gene Name Mapping.txt in the input folder.
                # I freely admit this was horribly designed though

                if strain == 'NONE':
                    name_dict[(gene_colloquial_name, species, None)] = mg1655_accession
                    reverse[mg1655_accession].append((gene_colloquial_name, species, None))
                else:
                    name_dict[(gene_colloquial_name, species, strain)] = mg1655_accession
                    reverse[mg1655_accession].append((gene_colloquial_name,
                                                      species,
                                                      constants.get_strain_converter(strain)))

            if len(exceptions) > 0:
                output_string = ''

                for (gene, strain) in exceptions:
                    output_string += str((gene, strain)) + '\n'

                raise AssertionError('Duplicate entries detected: %s' % output_string)

            return name_dict, reverse

        except:
            raise

    def revert(self, converted_value, throw_error=False):

        """
        
        Attempts to convert a standardized identifier back to the original IDs that match. It is not possible
        to recover the exact ID, but a set of possible original IDs will be returned.
        
        :param converted_value: str, identifier to revert
        :param throw_error: flag indicating whether or not to throw an error if an identifier is not found (true/false)
        :return: list of possible original (identifier, species) tuples correpsonding to converted_value
        """

        if converted_value in self.reversion:
            return self.reversion[converted_value]
        elif converted_value.upper() in self.reversion:
            return self.reversion[converted_value.upper()]
        elif not throw_error:
            return [(converted_value, 'unknown')]
        else:
            raise KeyError('Missing value from reversion dictionary: %s' % converted_value)

    def convert(self, name, species, strain=None, throw_error=False, remove=None):

        """
        
        Converts an identifier into its standardized form (if present).
        
        :param name: str, identifier
        :param species: str, associated species (or a filler string)
        :param strain: str, if identifier is associated with a specific species, None otherwise
        :param throw_error: boolean, throw error if name is missing 
        :param remove: str, strip off value of remove from an identifier, None if return identifiers as is
        :return: converted identifier
        """

        tag = ''

        if strain is not None:
            strain = strain.upper()
            strain = constants.get_strain_converter(strain).upper()

        if remove is not None and remove in name:
            name = name.replace(remove, '')
            tag = remove

        if name is None or species is None:
            return None

        if (name.upper(), species.upper(), strain) in self.standard:
            return self.standard[(name.upper(), species.upper(), strain)].upper() + tag
        if not throw_error:
            return name.upper() + tag

        raise KeyError('Name %s not found in standardization dictionary.' % name)

    def overwrite(self, standard_obj):

        """
        
        Overwrites entries in self.standard with those in standard_obj. This is used to ensure every identifier
        (as far as possible) has a standardized identifier by overwriting species specific entries with those
        that are more general. Anything remaining after that is the least specific possible identifier.
        
        :param standard_obj: standard object
        :return: None
        """

        for (name, species, strain) in standard_obj.standard:
            if (name, species, strain) in self.standard:
                self.standard[(name, species, strain)] = standard_obj.standard[(name, species, strain)]
