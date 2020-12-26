from resistome.utils import database_utils
from resistome.sql.validator import validate_mutation_data
import os
import psycopg2
import psycopg2.extras
from resistome import constants
from resistome.sql.sql_build_utils import prepare_sql_query, RESISTOME_SCHEMA, prepare_tuple_style_sql_query
from functools import lru_cache

__author__ = 'jdwinkler'

"""

This dict contains the column names for each table in the resistome.resistome postgres schema. It should be kept
consistent with the resistome_sql_schema.sql file in ../inputs . No type annotation is required as inserts are auto-
matically coerced to the type specified in the table definition with the exception of the json dicts used for 
encapsulating annotation data.

"""

annotation_format_rules = database_utils.parse_annotation_format_rules()
subobject_term_rules = database_utils.parse_input_terms()

position_cache = dict()
large_cache = dict()


def fetch_gene_position(cur, original_strain, converted_strain: str, gene_name: str):

    if (original_strain, converted_strain, gene_name) in position_cache:
        return position_cache[(original_strain, converted_strain, gene_name)]

    # get gene start from biocyc (#converted_strain#.genes) where converted_strain is
    # also the name of a biocyc table.
    if original_strain.lower() in constants.SPECIES_LIST or \
            constants.MAP_SPECIES_TO_REF.get(original_strain.lower(), 'Not in here!') in constants.SPECIES_LIST:
        cur.execute('select gene_locations.start, gene_locations.stop FROM genes '
                    'INNER JOIN gene_locations on gene_locations.gene_id = genes.gene_id '
                    'INNER JOIN strain on strain.strain_id = genes.strain_id '
                    'where upper(accession) = %s '
                    'AND strain.strain = lower(%s)',
                    (gene_name.upper(), converted_strain))
        gene_location = cur.fetchone()
    else:
        gene_location = None

    if gene_location is None:
        gene_start = 0
        gene_end = 0
    else:
        gene_start = int(gene_location['start'])
        gene_end = int(gene_location['stop'])

    position_cache[(original_strain, converted_strain, gene_name)] = (gene_start, gene_end, gene_location)

    return gene_start, gene_end, gene_location


def fetch_inbetween_genes(cursor, strain: str, accession_1: str, accession_2: str):
    """

    Fetches the genes that are in between gene1 and gene2 for the given strain. Results are inclusive of end genes
    specified if they exist. If the accession1 or accession2 are not found, then an array of [gene1, gene2] is returned.

    :param cursor:
    :param strain:
    :param accession_1:
    :param accession_2:
    :return:
    """

    if (strain, accession_1, accession_2) in large_cache:
        return large_cache[(strain, accession_1, accession_2)]

    if (strain, accession_2, accession_1) in large_cache:
        return large_cache[(strain, accession_2, accession_1)]

    cursor.execute('select accession from genes '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'INNER JOIN gene_locations on gene_locations.gene_id = genes.gene_id '
                   'WHERE strain.strain = %s '
                   'AND gene_locations.start >= (select min(gene_locations.start) from genes '
                   '                            INNER JOIN gene_locations ON gene_locations.gene_id = genes.gene_id '
                   '                            WHERE accession = ANY(%s)) '
                   'AND gene_locations.stop <= (select max(gene_locations.stop) from genes '
                   '                            INNER JOIN gene_locations ON gene_locations.gene_id = genes.gene_id '
                   '                            WHERE accession = ANY(%s)) '
                   'ORDER BY gene_locations.start ASC',
                   (strain, [accession_1, accession_2], [accession_1, accession_2]))

    output_genes = [x['accession'] for x in cursor.fetchall()]

    if accession_1 not in output_genes or accession_2 not in output_genes:
        # basically we could not locate on gene or the other. We should at least return the stated genes,
        # even if they are wrong.
        # note we do not check to make sure the provided boundaries are the actual boundaries.
        output = [accession_1, accession_2]
    else:
        output = output_genes

    large_cache[(strain, accession_1, accession_2)] = output

    return output


def insert_sql_get_id(schema, table, columns, id_name):
    """

    Appends 'returning ' + the name of the requested column. This is used in sql_generator to provide a way for
    a method to get the ID associated with an insert for subsequent related inserts.

    :param schema: name of the schema (str), usually resistome
    :param table: str, name of the table
    :param columns: iterable of strs, get from column dict
    :param id_name: name of column to return the value of post insert
    :return: formatted SQL query
    """

    sql = prepare_sql_query(schema=schema, table=table, columns=columns)

    sql += ' returning ' + id_name

    return sql


def paper_tag_table(cur, paper_id, tags):
    """

    Creates the paper tag table (paper id, 'tag' where tag is solvents_biofuels, metal_ions-a category of stress
    phenotype. One ID => many tags.

    :param cur: database cursor to the resistome db
    :param paper_id: unique ID (generated after insert into paper table) for the study
    :param tags: iterable of str with the resistome categorization of the study
    :return: None
    """

    sql = prepare_sql_query(schema='resistome', table='paper_tags', columns=RESISTOME_SCHEMA['paper_tags'])

    for tag in tags:
        cur.execute(sql, (paper_id, tag.lower()))


def error_check_annotation(mutation_type, annotation_data):
    """

    Checks annotation formats and content for validity.

    Ensures bases, residues are valid entries (defined AA/nucleotide bases in database_contents
    Ensures format is consistent for classes of entries

    Returns true if a given mutation type, annotation_data pairing are valid, False otherwise.

    :param mutation_type:
    :param annotation_data:
    :return:

    """

    passes_error_check = True
    value = annotation_data[mutation_type]

    allowed_position_types = {'absolute_inferred', 'absolute', 'gene_unknown'}

    try:
        if mutation_type == 'compartmentalization':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'aa_snps':
            for [position, original_aa, mutated_aa] in value:
                if original_aa not in constants.AMINO_ACIDS or mutated_aa not in constants.AMINO_ACIDS:
                    passes_error_check = False
                if position is None or not isinstance(position, int):
                    passes_error_check = False
        elif mutation_type == 'nuc_snps':
            for x in value:
                (position, original_base, mutated_base) = (x[0], x[1], x[2])
                if not isinstance(position, int) and position is not None:
                    passes_error_check = False
                if original_base not in constants.NUCLEOTIDE_BASES or mutated_base not in constants.NUCLEOTIDE_BASES:
                    passes_error_check = False
        elif mutation_type == 'indel':
            for x in value:
                (position, length, position_type) = (x[0], x[1], x[2])
                length = int(length)
                if not isinstance(position, int) and position is not None:
                    passes_error_check = False
                if not isinstance(length, int):
                    passes_error_check = False
                if position_type not in allowed_position_types:
                    passes_error_check = False
                if isinstance(length, int) and length == 0:
                    passes_error_check = False
        elif mutation_type == 'intergenic':
            first_gene = value[0]
            second_gene = value[1]
            passes_error_check = True if isinstance(first_gene, str) and isinstance(second_gene, str) else False
        elif mutation_type == 'is_insertion':
            (insertion_element, descriptor) = (value[0], value[1])
            if not isinstance(insertion_element, str):
                passes_error_check = False
        elif mutation_type == 'amplified':
            passes_error_check = True if isinstance(value, int) else False
        elif mutation_type == 'antisense':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'scaffold_binder':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'truncated':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'codonoptimized':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'con':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'large_amplification':
            (first_gene, last_gene, degree_of_amplification) = value
            passes_error_check = True if isinstance(first_gene, str) and \
                                         isinstance(last_gene, str) and \
                                         isinstance(degree_of_amplification, int) else False
        elif mutation_type == 'large_deletion' or mutation_type == 'large_inversion':
            (first_gene, last_gene) = value
            passes_error_check = True if isinstance(first_gene, str) and \
                                         isinstance(last_gene, str) \
                else False
        elif mutation_type == 'del':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'frameshift':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'protein_fusion':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'duplication':
            for x in value:
                (start, stop, length, position_type) = (x[0], x[1], x[2], x[3])
                if not isinstance(start, int):
                    passes_error_check = False
                if not isinstance(stop, int):
                    passes_error_check = False
                if not isinstance(length, int):
                    passes_error_check = False
                if position_type not in allowed_position_types:
                    passes_error_check = False
                if isinstance(length, int) and length == 0:
                    passes_error_check = False
        elif mutation_type == 'mutated':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'oe':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'plasmid':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'promoter_switch':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'rbs_tuned':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'less_nadh_inhibition':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'regulated':
            for x in value:
                passes_error_check = True if isinstance(x, str) else False
                if not passes_error_check:
                    break
        elif mutation_type == 'rep':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'scaffold_bindee':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'mrna_secondarystructure':
            passes_error_check = True if isinstance(value, str) or value is None else False
        elif mutation_type == 'sensor':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'integrated':
            passes_error_check = True if isinstance(value, str) else False
        elif mutation_type == 'terminated':
            passes_error_check = True if isinstance(value, str) else False
        else:
            raise AssertionError('Unknown mutation type: %s' % mutation_type)

    except:
        passes_error_check = False

    return passes_error_check


def error_check_mutation_entry(entry_data, mutation_types):
    """
    Checks to ensure that all required fields are defined in mutation data, and that they conform to the expected
    format specified in Input Terms_typed.txt.

    :param entry_data: tuple of data constructed by insert_mutation_data function
    :param mutation_types: list of mutations encoded in the mutation object
    :return:
    """

    (paper_id, unique_mutant_id, gene_name, species, converted_strain, effects, original) = entry_data

    passes_error_check = True

    if gene_name is None:
        passes_error_check = False
    if ',' in gene_name or '/' in gene_name:
        passes_error_check = False
    if 'Escherichia coli'.upper() not in species.upper():
        passes_error_check = False
    if converted_strain not in constants.SPECIES_LIST:
        passes_error_check = False

    for change in mutation_types:
        if change not in annotation_format_rules:
            passes_error_check = False

    return passes_error_check


def error_check_expression_entry(entry_data):
    """

    Checks to ensure that all required fields are defined in expression data, and that they conform to the expected
    format specified in Input Terms_typed.txt.

    :param entry_data:
    :return:
    """

    (paper_id, unique_mutant_id, study_id, gene_name, species, converted_strain, ge_change, p_value,
     fold_change) = entry_data

    passes_error_check = True

    if gene_name is None:
        passes_error_check = False
    if 'Escherichia coli'.upper() not in species.upper():
        passes_error_check = False
    if converted_strain not in constants.SPECIES_LIST:
        passes_error_check = False
    if ge_change not in {'+', '-'}:
        passes_error_check = False
    if p_value is not None and not isinstance(p_value, float):
        passes_error_check = False
    if fold_change is not None and not isinstance(fold_change, float):
        passes_error_check = False

    return passes_error_check


def error_check_expression_study(entry_data):
    """

    Checks to ensure that all required fields are defined in expression study data, and that they conform to the expected
    format specified in Input Terms_typed.txt.

    :param entry_data:
    :return:
    """

    (paper_id,
     unique_mutant_id,
     accession,
     ge_method,
     stat_test,
     exposure,
     growth_phase,
     stressor,
     units) = entry_data

    passes_error_check = True

    if growth_phase is not None and growth_phase not in constants.GROWTH_PHASES:
        passes_error_check = False
    if ge_method is not None and ge_method not in constants.EXPRESSION_METHOD:
        passes_error_check = False

    return passes_error_check


def error_check_mutant_entry(entry_data, methods):
    """

    Checks to ensure that all required fields are defined in mutant data, and that they conform to the expected
    format specified in Input Terms_typed.txt.

    :param entry_data:
    :return:
    """

    (paper_id,
     mutant_name,
     species,
     strain,
     converted_strain,
     oxygen,
     medium,
     carbon,
     supplements,
     ph,
     cvessel,
     vvolume,
     cvolume,
     temperature,
     rotation,
     fold_imp,
     init_fit,
     final_fit,
     fit_unit,
     comments) = entry_data

    passes_error_check = True

    if 'Escherichia coli'.upper() not in species.upper():
        passes_error_check = False
    if converted_strain not in constants.SPECIES_LIST:
        passes_error_check = False
    if oxygen is not None and oxygen not in constants.AIR:
        passes_error_check = False
    if carbon is None:
        passes_error_check = False
    if ph is not None and not isinstance(ph, float):
        passes_error_check = False
    if rotation is not None and not isinstance(rotation, float):
        passes_error_check = False
    if cvolume is not None and not isinstance(cvolume, float):
        passes_error_check = False
    if vvolume is not None and not isinstance(vvolume, float):
        passes_error_check = False
    if temperature is not None and not isinstance(temperature, float):
        passes_error_check = False
    if cvessel is not None and cvessel not in constants.CULTURE_VESSELS:
        passes_error_check = False

    for method in methods:
        if method not in constants.ENGINEERING_METHODS:
            passes_error_check = False

    return passes_error_check


def error_check_paper_entry(entry_data):
    """

    Checks to ensure that all required fields are defined in paper data, and that they conform to the expected
    format specified in Input Terms_typed.txt.

    Note: the duplicate DOI check is located in the ContainerClass method in object_representation.py.

    :param entry_data:
    :return:
    """

    (title, doi, year, group, journal, method, difficulty, reason, reference_genome, total_designs,
     comments) = entry_data

    passes_error_check = True

    if title is None or doi is None:
        passes_error_check = False
    if doi is not None and len(doi) == 0:
        passes_error_check = False
    if journal is None:
        passes_error_check = False
    if group is None:
        passes_error_check = False
    if not isinstance(year, int):
        passes_error_check = False
    if journal is None:
        passes_error_check = False
    if method is None or not set(method).issubset(constants.EXPERIMENTAL_DESIGNS):
        passes_error_check = False
    if not isinstance(difficulty, int) or (difficulty < 0 or difficulty > 5):
        passes_error_check = False
    if reason is None or not set(reason).issubset(constants.SCORE_REASONS):
        passes_error_check = False

    return passes_error_check


def insert_mutant_phenotype(cur, unique_mutant_id, phenotype, type_of_phenotype, r_level, r_units, tag_dict, root_dict):
    """

    Prepares inserts for the phenotype table in the resistome. This table links the mutant entry to data about the
    phenotype(s) in question, including the resistance level and units.

    Tag dict and root dict are defined in the ../inputs/ directory in the Stress Classification.txt and
    Compound Information.txt. It is assumed that every phenotype/compound encountered in the database is listed
    in both files, otherwise an error is thrown.

    :param cur: database cursor for the resistome
    :param unique_mutant_id: unique mutant id after insert into the mutant table
    :param phenotype: str, name of phenotype
    :param type_of_phenotype: str, 'R' for resistant or 'S' for sensitive
    :param r_level: str, some type of representation about the resistance level (often numeric)
    :param r_units: str, units of r_level
    :param tag_dict: dict (str: str) converts phenotype to its general class
    :param root_dict: dict (str, tuple) describes highest level class of the phenotype
    :return: None
    """

    sql = prepare_sql_query(schema='resistome', table='phenotypes',
                            columns=RESISTOME_SCHEMA['phenotypes'])

    error_phenotypes = []

    if phenotype not in tag_dict:
        error_phenotypes.append('Missing phenotype name: %s' % (phenotype,))

    if phenotype is None or len(phenotype) == 0:
        error_phenotypes.append('Missing phenotype?')

    cur.execute(sql, (unique_mutant_id,
                      phenotype.lower(),
                      tag_dict.get(phenotype, phenotype).lower(),
                      type_of_phenotype,
                      root_dict.get(phenotype, 'X')[0],
                      r_level,
                      r_units))

    return error_phenotypes


def insert_annotation(cur, unique_gene_id, annotation, json_dict):
    """

    Inserts a json version of the annotations for each mutation into the annotation table. Json is actually a terrible
    type for doing this, but given that annotations take a wide variety of forms and the mutation data have to be
    handled on a case by case basis anyway, this construction seems to work well enough for now.

    Unfortunately the Resistome is gene-centric, so this representation gets a bit awkward in case large non-coding
    regions are mutated.

    :param cur: database cursor for the resistome
    :param unique_gene_id: unique id for the associated mutation entry
    :param annotation: str, resistome type of mutation
    :param json_dict: dict of (str, mutation type: (data)) to contain the annotation
    :return: None
    """

    error_check = error_check_annotation(mutation_type=annotation,
                                         annotation_data=json_dict)

    return error_check, (unique_gene_id, annotation.lower(), psycopg2.extras.Json(json_dict))


def insert_mutation_data(cur, year, paper_id, unique_mutant_id, mutation_obj, species, converted_strain,
                         original_strain):
    """

    Inserts mutation data from Resistome records into the mutations table.

    :param cur: database cursor for the resistome
    :param paper_id: unique id of the associated study
    :param unique_mutant_id: unique id of the associated mutant
    :param mutation_obj: a mutation object created by MetEngDatabase when parsing the corresponding .txt file
    :param species: str, species as written in the original record
    :param converted_strain: str, converted species name (REL606, MG1655, W)
    :return: None
    """

    original = mutation_obj.is_original
    effects = mutation_obj.effects
    changes = [x.lower() for x in mutation_obj.changes]

    output_errors = []

    # genes that will be inserted into the table
    # can be quite long for large_(amplification/deletion/inversion) mutations
    genes_to_add = list(mutation_obj.modified_genes)

    if len(changes) == 0:
        output_errors.append('Missing mutation types (GeneMutation field)?')

    large_mutation_boundaries = None

    if 'large_amplification' in changes:
        changes = ['large_amplification']
        try:
            # converts the mutation into gene1, geneN, amount of amplification
            mutation_obj.annotation_backing['large_amplification'] = (mutation_obj.annotation('large_amplification')[0],
                                                                      mutation_obj.annotation('large_amplification')[-2],
                                                                      int(mutation_obj.annotation(
                                                                          'large_amplification')[-1]))
            large_mutation_boundaries = (mutation_obj.annotation('large_amplification')[0],
                                         mutation_obj.annotation('large_amplification')[-2])
        except (ValueError, TypeError):
            # this sometimes occurs if the value in the file is actually a float
            # though a sane language would probably just truncate the .0...
            mutation_obj.annotation_backing['large_amplification'] = (mutation_obj.annotation('large_amplification')[0],
                                                                      mutation_obj.annotation('large_amplification')[-2],
                                                                      int(float(mutation_obj.annotation(
                                                                          'large_amplification')[-1])))
            large_mutation_boundaries = (mutation_obj.annotation('large_amplification')[0],
                                         mutation_obj.annotation('large_amplification')[-2])

    elif 'large_deletion' in changes:

        changes = ['large_deletion']
        # converts the mutation into gene1, geneN
        mutation_obj.annotation_backing['large_deletion'] = (mutation_obj.annotation_backing['large_deletion'][0],
                                                             mutation_obj.annotation_backing['large_deletion'][-1])
        large_mutation_boundaries = (mutation_obj.annotation('large_deletion')[0],
                                     mutation_obj.annotation('large_deletion')[-1])

    elif 'large_inversion' in changes:
        # converts the mutation into gene1, geneN
        mutation_obj.annotation_backing['large_inversion'] = (mutation_obj.annotation_backing['large_inversion'][0],
                                                              mutation_obj.annotation_backing['large_inversion'][-1])
        large_mutation_boundaries = (mutation_obj.annotation('large_inversion')[0],
                                     mutation_obj.annotation('large_inversion')[-1])

    if large_mutation_boundaries is not None:
        genes_to_add = fetch_inbetween_genes(cursor=cur, strain=converted_strain,
                                             accession_1=large_mutation_boundaries[0],
                                             accession_2=large_mutation_boundaries[1])

    upload_tuples = []

    for gene_name in genes_to_add:

        data_tuple = (paper_id,
                      unique_mutant_id,
                      gene_name.upper().replace('\'', ''),
                      species,
                      converted_strain.lower(),
                      effects,
                      original)

        valid_mutation_object = error_check_mutation_entry(data_tuple, changes)

        if not valid_mutation_object:
            output_errors.append('Error in mutation object: %s, %s' % (str(data_tuple), str(changes)))

        # get unique id to associate the annotation with the proper mutation
        sql = insert_sql_get_id('resistome', 'mutations', RESISTOME_SCHEMA['mutations'], 'gene_id')
        cur.execute(sql, data_tuple)
        gene_id = cur.fetchone()['gene_id']

        # get gene start from biocyc (#converted_strain#.genes) where converted_strain is
        # also the name of a biocyc table.
        gene_start, gene_end, gene_location = fetch_gene_position(cur, original_strain, converted_strain,
                                                                  gene_name)

        proper_location_tuple = constants.get_proper_start_stop(gene_name.upper(),
                                                                year,
                                                                query_species=converted_strain.lower())

        if proper_location_tuple is not None:
            (ref_gene_start, ref_gene_end) = proper_location_tuple
        else:
            ref_gene_start = 0
            ref_gene_end = 0

        for change in changes:

            annotation = mutation_obj.annotation(change)

            if change == 'nuc_snps':
                # all locations are absolute by default
                # todo: add in genome reference correction based on year (2013)
                annotation = [x.replace('|absolute', '') for x in annotation]

                annotation_tuples = []

                for entry in annotation:

                    original_value = entry[0]
                    altered_value = entry[-1]
                    try:
                        location = int(entry[1:-1])
                    except:
                        if entry[1:-1] != '?':
                            raise
                        else:
                            location = None

                    # Update DEC2020: after thinking it over when adding new data, it does not make sense to correct
                    # for relative/absolute locations.

                    if gene_start != 0 and gene_end != 0 and isinstance(location, int):
                        relative_location = constants.convert_to_relative_position(gene_start,
                                                                                   gene_end,
                                                                                   location)
                        location = constants.convert_to_absolute_position(gene_start,
                                                                          gene_end,
                                                                          relative_location)

                    annotation_tuples.append((location, original_value, altered_value))

                annotation = annotation_tuples

            if change == 'aa_snps':

                annotation_tuples = []

                for entry in annotation:

                    original_value = entry[0].upper()
                    altered_value = entry[-1].upper()
                    location = entry[1:-1]

                    if 'STOP' in location:
                        raise AssertionError('Found STOP in codon, should replace with *!')

                    try:
                        # fix zero indexing
                        location = int(location) - 1
                        pass
                    except (ValueError, TypeError):
                        # if location is not actually given (very rare)
                        location = None
                    annotation_tuples.append((location, original_value, altered_value))
                annotation = annotation_tuples

            if change == 'indel':
                output_tuples = []
                for (location, size, type_of_reference) in annotation:
                    if type_of_reference == 'relative' or type_of_reference is None:
                        if gene_location is None:
                            if location == 'coding':
                                location = None
                            output_tuples.append((location, size, 'gene_unknown'))
                            continue

                        if ref_gene_start != 0 and ref_gene_end != 0:
                            # updated reference
                            actual_location = constants.convert_to_absolute_position(gene_start, gene_end, location)
                            location_type = 'absolute'
                        elif gene_start != 0 and gene_end != 0 and location != 'coding':
                            # new reference, actual location provided
                            actual_location = constants.convert_to_absolute_position(gene_start, gene_end, location)
                            location_type = 'absolute'
                        elif gene_start != 0 and gene_end != 0 and location == 'coding':
                            # new reference, actual location not provided
                            actual_location = int((float(gene_start) + float(gene_end)) / 2.0)
                            location_type = 'absolute_inferred'
                        else:
                            # should be no unhandled cases
                            raise AssertionError('Should not exist.')
                        output_tuples.append((actual_location, size, location_type))

                    else:
                        if location == 'coding':
                            location = None
                        output_tuples.append((location, size, type_of_reference))

                annotation = output_tuples

            if change == 'duplication':

                output_tuples = []

                (location, size, type_of_reference) = annotation
                location_tokens = location.split('-')
                start = int(location_tokens[0])
                stop = int(location_tokens[1])
                size = int(size)

                if type_of_reference == 'relative' or type_of_reference is None:

                    if gene_location is None:
                        output_tuples.append((location, size, 'gene_unknown'))
                        continue

                    if gene_end != 0 and gene_start != 0:
                        actual_start_location = constants.convert_to_absolute_position(gene_start, gene_end, start)
                        actual_stop_location = constants.convert_to_absolute_position(gene_start, gene_end, stop)
                        location_type = 'absolute'
                        output_tuples.append((actual_start_location, actual_stop_location, size, location_type))
                    else:
                        output_tuples.append((start, stop, size, type_of_reference))
                else:
                    output_tuples.append((start, stop, size, type_of_reference))

                annotation = output_tuples

            if change == 'oe':
                annotation = annotation[0]

            if change == 'plasmid':
                annotation = annotation[0]

            wrapper_dict = {change: annotation}
            valid, tuple_to_upload = insert_annotation(cur, gene_id, change, wrapper_dict)
            upload_tuples.append(tuple_to_upload)

            if not valid:
                output_errors.append(('Error detected in %s: %s (%s)' % (change, str(wrapper_dict[change]),
                                                                         gene_name)))

    sql = prepare_tuple_style_sql_query(schema='resistome',
                                        table='annotations',
                                        columns=RESISTOME_SCHEMA['annotations'])

    psycopg2.extras.execute_values(cur=cur, sql=sql, argslist=upload_tuples, page_size=2000)

    return output_errors


def insert_expression_data(cur, paper_id, unique_mutant_id, study_id, expression_obj, species, converted_strain):
    """

    Inserts expression data into the expressions table. Conceptually very similar to the mutations table,
    but without a companion 'annotations' table.

    :param cur: database cursor for the resistome
    :param paper_id: unique id for the study
    :param unique_mutant_id: unique id for the mutant
    :param study_id: unique id for the expression study
    :param expression_obj: an object describing the fields associated with annotated resistome expression data
    :param species: original species name specified in Resistome record
    :param converted_strain: converted strain name (mg1655, rel606, w)
    :return:
    """

    errors = []

    gene_name = expression_obj.name
    try:
        ge_change = '+' if expression_obj.gechange.lower() == 'overexpressed' else '-'
    except:
        # this occurs if gechange is missing
        if expression_obj.foldchange > 1.0:
            ge_change = '+'
        else:
            ge_change = '-'
    p_value = expression_obj.pvalue
    fold_change = expression_obj.foldchange

    data_tuple = (paper_id,
                  unique_mutant_id,
                  study_id,
                  gene_name,
                  species,
                  converted_strain,
                  ge_change,
                  p_value,
                  fold_change)

    valid_expression_object = error_check_expression_entry(data_tuple)

    if not valid_expression_object:
        errors.append('Error in expression object data: %s' % str(data_tuple))

    sql = prepare_sql_query(schema='resistome', table='expressions', columns=RESISTOME_SCHEMA['expressions'])
    cur.execute(sql, data_tuple)

    return errors


def insert_mutant(cur, year, paper_id, mutant_obj, tag_info, root_classes):
    """

    Inserts a mutant record into the mutants table in the resistome. This method does the lion's share of parsing,
    as it drives inserts into the phenotype, methods, mutations, and expression_studies/expression tables
    since they are all based on connections to the unique mutant genotype being inserted here.

    :param cur: database cursor for the resistome
    :param paper_id: unique id associated with the study
    :param mutant_obj: object created to represent a given mutant entry by MetEngDatabase
    :param tag_info: dict (str : str) to annotate specific phenotypes with their general class
    :param root_classes: dict (str : tuple) higher level hierarchy to describe phenotypes
    :return: None
    """

    errors = []

    mutant_name = mutant_obj.name
    species = mutant_obj.species
    strain = mutant_obj.strain.upper()

    if strain == 'K12':
        strain = 'K-12'
    if strain == 'K':
        strain = 'K-12'

    # culture vessel, options listed in Term Usage.txt under ../inputs/
    cvessel = mutant_obj.culture_vessel

    # oxygenation
    oxygen = mutant_obj.oxygenation
    (medium, supplements) = mutant_obj.medium
    carbon = mutant_obj.substrates
    ph = mutant_obj.pH
    cvolume = mutant_obj.culture_volume
    vvolume = mutant_obj.liquid_volume
    temperature = mutant_obj.temperature
    rotation = mutant_obj.rpm
    fold_imp = mutant_obj.fold_improvement
    init_fit = mutant_obj.initial_fitness
    final_fit = mutant_obj.final_fitness
    fit_unit = mutant_obj.fitness_unit

    # comments is a free-form text field not used for anything currently
    comments = mutant_obj.comments

    methods = mutant_obj.methods

    resist_level = mutant_obj.resistance_level
    resist_units = mutant_obj.resistance_units

    sensitive_phenotypes = mutant_obj.sensitive_phenotypes
    resistant_phenotypes = mutant_obj.resistant_phenotypes

    converted_strain = constants.get_strain_converter(strain)

    mutant_tuple = (paper_id,
                    mutant_name,
                    species,
                    strain,
                    converted_strain,
                    oxygen,
                    medium,
                    carbon,
                    supplements,
                    ph,
                    cvessel,
                    vvolume,
                    cvolume,
                    temperature,
                    rotation,
                    fold_imp,
                    init_fit,
                    final_fit,
                    fit_unit,
                    comments)

    valid_mutant_object = error_check_mutant_entry(mutant_tuple, methods)

    if not valid_mutant_object:
        errors.append('Error detected in high-level mutant object: %s, %s' % (str(mutant_tuple), str(methods)))

    sql = insert_sql_get_id('resistome', 'mutants', RESISTOME_SCHEMA['mutants'], 'mutant_id')
    cur.execute(sql, mutant_tuple)
    unique_mutant_id = cur.fetchone()['mutant_id']

    study_id = None

    if len(mutant_obj.expression) > 0:
        expression_obj = mutant_obj.expression[0]

        stat_test = expression_obj.stattest
        exposure = expression_obj.exposure
        growth_phase = expression_obj.growthphase
        ge_method = expression_obj.gemethod
        stressor = expression_obj.stress_amount
        units = expression_obj.stress_units
        accession = expression_obj.accession

        if growth_phase is not None:
            growth_phase = growth_phase.lower()

        if ge_method is not None:
            ge_method = ge_method.lower()

        high_level_descriptor = (paper_id,
                                 unique_mutant_id,
                                 accession,
                                 ge_method,
                                 stat_test,
                                 exposure,
                                 growth_phase,
                                 stressor,
                                 units)

        valid_expression_study = error_check_expression_study(high_level_descriptor)

        if not valid_expression_study:
            errors.append('Error in expression study: %s' % str(high_level_descriptor))

        sql = insert_sql_get_id('resistome', 'expression_studies', RESISTOME_SCHEMA['expression_studies'], 'study_id')
        cur.execute(sql, high_level_descriptor)
        study_id = cur.fetchone()['study_id']

    for method in methods:
        sql = prepare_sql_query(schema='resistome', table='mutant_methods',
                                columns=RESISTOME_SCHEMA['mutant_methods'])
        cur.execute(sql, (unique_mutant_id, method))

    for phenotype in sensitive_phenotypes:
        errors.extend(insert_mutant_phenotype(cur, unique_mutant_id, phenotype, 'S',
                                              resist_level, resist_units, tag_info, root_classes))

    for phenotype in resistant_phenotypes:
        errors.extend(insert_mutant_phenotype(cur, unique_mutant_id, phenotype, 'R',
                                              resist_level, resist_units, tag_info, root_classes))

    affected_genes = []

    for mutation in mutant_obj.mutations:
        affected_genes.extend(mutation.modified_genes)
        errors.extend(insert_mutation_data(cur, year, paper_id, unique_mutant_id, mutation, species, converted_strain,
                                           strain))

    for expression in mutant_obj.expression:
        errors.extend(insert_expression_data(cur, paper_id, unique_mutant_id, study_id, expression, species,
                                             converted_strain))

    return errors


def build_gene_standardization_table(cur, papers, species_std, unified_std):
    """

    Builds the gene standardization table (essentially tuples of E. coli species, accession, and the mg1655 equivalent
    if there is one).

    :param cur: database cursor for the resistome
    :param papers: collection of paper objects (paper with [mutant : [mutations]]
    :return: None
    """

    conversion_table_sql = prepare_tuple_style_sql_query(schema='resistome',
                                                         table='gene_standardization',
                                                         columns=RESISTOME_SCHEMA['gene_standardization'])
    metadata_table_sql = prepare_tuple_style_sql_query(schema='resistome', table='gene_metadata',
                                                       columns=RESISTOME_SCHEMA['gene_metadata'])

    added_already = set()

    in_standard = set(species_std.standard.values())

    from collections import defaultdict
    multiple_conversion = defaultdict(list)

    metadata_tuples = []
    conversion_tuples = []

    # extract standardization data from 'Species Specific Gene Names.txt' in ../inputs/
    for (gene, species, strain) in species_std.standard:
        standardized_name = species_std.convert(gene, species, strain=strain)
        mg1655_name = unified_std.convert(gene, species, strain='mg1655')

        multiple_conversion[(gene, strain)].append((standardized_name, mg1655_name))

        if (strain, standardized_name) not in added_already:
            conversion_tuples.append((strain, standardized_name, mg1655_name))

        metadata_tuples.append((gene, species, strain, standardized_name, mg1655_name))
        added_already.add((strain, standardized_name))

    for (gene, strain) in multiple_conversion:

        if len(multiple_conversion[(gene, strain)]) > 1:
            print(gene, strain, multiple_conversion[(gene, strain)])

    for paper in papers:
        for mutant in paper.mutants:
            strain = constants.get_strain_converter(mutant.strain)
            for mutation in mutant.mutations:
                for gene in mutation.modified_genes:
                    standardized_name = species_std.convert(gene, mutant.species, strain=strain)
                    # adds anything that isn't covered by the standardization file (assumes it is a weirdly named
                    # gene-not very common, but often enough to handle
                    if gene.upper() == standardized_name.upper() and gene.upper() not in in_standard and (
                    strain.upper(), gene) not in added_already:
                        metadata_tuples.append((gene, mutant.species, strain, gene, gene))
                        conversion_tuples.append((strain.upper(), gene, gene))
                        added_already.add((strain.upper(), gene))

    psycopg2.extras.execute_values(cur, metadata_table_sql, argslist=metadata_tuples, page_size=2000)
    psycopg2.extras.execute_values(cur, conversion_table_sql, argslist=conversion_tuples, page_size=2000)


def build_phenotype_standardization_table(cur):
    """

    Builds the phenotype standardization table compiling all phenotypes in the resistome and their associated
    categories/classifications.

    :param cur: database cursor for the resistome
    :return: None
    """

    phenotype_std = database_utils.get_phenotype_standard()

    tag_dict, ontology = database_utils.load_tolerance_ontology()

    pk = dict()
    tuples_to_upload = []

    tags_to_consider = set()
    for (tag, _, _) in phenotype_std.standard:
        tag = tag.lower()
        tags_to_consider.add(tag.replace('_sensitive', ''))

    for tag in tag_dict:
        tag = tag.lower()
        tags_to_consider.add(tag.replace('_sensitive', ''))

    for tag in tags_to_consider:

        standardized_tag = ''

        try:
            # try to convert the phenotype to something more standardized
            standardized_tag = phenotype_std.convert(tag, 'None')

            # if the above fails, tries the original tag (might already by the standard form)
            tag_category = tag_dict.get(standardized_tag, tag_dict.get(tag, None))

            # root_class is something like chemical or environmental stress, specific classes as you'd expect
            (root_class, specific_classes) = ontology.get(standardized_tag, ontology.get(tag, None))
        except:
            raise ValueError('Missing stress ontology tag: %s, %s' % (standardized_tag, tag))

        if standardized_tag in pk:
            tuples_to_upload.append((tag.lower(), pk[standardized_tag]))
            continue

        sql = prepare_sql_query(schema='resistome',
                                table='phenotype_standardization',
                                columns=RESISTOME_SCHEMA['phenotype_standardization']) + ' RETURNING phenotype_id'
        cur.execute(sql, (standardized_tag.lower(),
                          tag_category.lower(),
                          root_class.lower(),
                          list(specific_classes)))

        phenotype_id = cur.fetchone()['phenotype_id']
        tuples_to_upload.append((tag.lower(), phenotype_id))
        pk[standardized_tag] = phenotype_id

    psycopg2.extras.execute_values(cur,
                                   sql='INSERT INTO resistome.phenotype_mapping (original_name, standardized_id) VALUES %s',
                                   argslist=tuples_to_upload)


def build_term_explanation_table(cur, filepath: str):
    """

    Loads the term data originally used for the resistome local data entry tool into a table. This is mainly for
    user reference and is not used anywhere.

    :param cur: database cursor for the resistome
    :return: None
    """

    with open(filepath, 'r') as fhandle:
        for line in fhandle:
            tokens = line.strip().split('\t')
            position = tokens[0]
            internal_name = tokens[1]

            # public explanation is the converted text display to the user, so often it is less of an explanation
            # and more of an interrogatory statement
            public_explanation = tokens[2]

            sql = prepare_sql_query(schema='resistome', table='term_explanation',
                                    columns=RESISTOME_SCHEMA['term_explanation'])
            cur.execute(sql, (position, internal_name, public_explanation))


def build_abbreviation_tables(cur):
    """

    Centralizes the documents describing various forms of abbreviations for phenotypes, experimental methods,
    journals, mutations, and general phenotype categories.

    :param cur: database cursor for the resistome
    :return: None
    """

    files = [('category_abbrevs.txt', 'categories'),
             ('journal_abbrevs.txt', 'journal'),
             ('method_abbrevs.txt', 'methods'),
             ('mutation_abbrevs.txt', 'mutations'),
             ('tolerance_abbrevs.txt', 'phenotype')]

    for (file_x, table_type) in files:

        with open(os.path.join(constants.INPUT_DIR, 'settings', file_x)) as f:

            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                entry = tokens[0].upper()
                conversion = tokens[1].upper()

                sql = prepare_sql_query(schema='resistome', table='abbreviations',
                                        columns=RESISTOME_SCHEMA['abbreviations'])
                cur.execute(sql, (entry.upper(), table_type, conversion.upper()))


def build_go_metabolite_tables(cur):
    """

    Creates two tables containing gene accession to either metabolite (based on Zampieri et al 2017, MSB, Metabolic
    Constraints on Antibiotic Resistance) or gene ontology annotations per gene from Biocyc.

    :param cur: database cursor for the resistome
    :return: None
    """

    with open(os.path.join(constants.INPUT_DIR, 'metabolomics', 'mg1655_accession_to_metabolite.txt'), 'r') as f:

        for line in f:
            tokens = line.strip().split('\t')
            accession = tokens[0]
            metabolite = tokens[1]

            # note that the Zampieri study only has MG1655 data, and it is probably a bad idea to extrapolate
            # to other subspecies without explicit user acknowledgement of that fact
            sql = prepare_sql_query(schema='resistome',
                                    table='metabolomics',
                                    columns=RESISTOME_SCHEMA['metabolomics'])
            cur.execute(sql, (accession, metabolite))

    with open(os.path.join(constants.INPUT_DIR, 'ontologies', 'accession_to_go_terms.txt'), 'r') as f:

        for line in f:
            tokens = line.strip().split('\t')
            accession = tokens[0]
            metabolite = tokens[1]
            sql = prepare_sql_query(schema='resistome',
                                    table='gene_ontology',
                                    columns=RESISTOME_SCHEMA['gene_ontology'])
            cur.execute(sql, (accession, metabolite))


def main(cur, add_resistome_go_metabolite_tables=False):

    from collections import defaultdict

    species_std = database_utils.get_species_standard()
    unified_std = database_utils.get_unified_standard()

    tags_information, ontology = database_utils.load_tolerance_ontology()

    # other databases currently only includes Carol Gross' large phenotyping effort for the Keio collection
    # NOTE THIS IS WHERE YOU CONTROL THE INPUTS
    ### ###

    ### ###
    # Can exclude files if you want!
    standardized_papers = database_utils.get_database(folders=('2018', '2017', '2016', '2015', 'other_databases',
                                                               '2019', '2020'))
    print('Inserting Resistome Data into SQL database')

    build_gene_standardization_table(cur, standardized_papers, species_std=species_std, unified_std=unified_std)

    # many tables link against here to ensure that they are being consistent with Resistome standards.
    build_term_explanation_table(cur, filepath=os.path.join(constants.INPUT_DIR, 'settings', 'FK_InternalFields.txt'))

    build_phenotype_standardization_table(cur)

    build_abbreviation_tables(cur)

    # if add_resistome_go_metabolite_tables:
    #     build_go_metabolite_tables(cur)

    error_dict = defaultdict(set)

    doi_duplicate_tracker = set()

    for k, paper in enumerate(standardized_papers):

        print('(%04i/%04i) Working on %s [Tags=%s]' % (k, len(standardized_papers), paper.title,
                                                       ', '.join(paper.categories)))

        title = paper.title
        doi = paper.doi
        year = paper.year
        group = paper.rgroup
        journal = paper.journal
        difficulty = paper.difficulty
        reason = paper.difficulty_reason
        method = paper.design_method
        tags = paper.categories
        total_designs = paper.total_designs
        comments = paper.comments

        # this check is in MetEngDatabase as well
        if doi in doi_duplicate_tracker:
            error_dict[title].add('Duplicate DOI detected: %s' % doi)

        if difficulty is None:
            difficulty = 1

        # human = none based on Terms Usage.txt
        if 'human' in method:
            method = set(method)
            method.remove('human')
            method.add('none')
            method = list(method)

        reference_genome = paper.reference_genome()

        if total_designs is None:
            total_designs = len(paper.mutants)

        paper_tuple = (title, doi, year, group, journal, method, difficulty, reason, reference_genome, total_designs,
                       comments)

        valid_paper_object = error_check_paper_entry(paper_tuple)

        if not valid_paper_object:
            error_dict[title].add('Error in paper level object: %s' % str(paper_tuple))

        sql = insert_sql_get_id('resistome',
                                'papers',
                                RESISTOME_SCHEMA['papers'],
                                'paper_id')

        cur.execute(sql, paper_tuple)
        inserted_id = cur.fetchone()['paper_id']

        paper_tag_table(cur, inserted_id, tags)

        for mutant in paper.mutants:

            paper_errors = insert_mutant(cur, year, inserted_id, mutant, tags_information, ontology)
            if len(paper_errors) > 0:
                error_dict[title].update(paper_errors)

    if len(error_dict.values()) > 0:

        for title in error_dict:
            if len(error_dict[title]) > 0:
                print('Paper title: %s' % title)
                for e in error_dict[title]:
                    print(e)

        raise AssertionError('Errors detected in data corpus; aborting insert transaction.')

    print('Beginning data validation.')

    invalid_annotation_entries, count = validate_mutation_data(cur)

    print('Got (in)validated set of annotations. Altering annotation table valid, reason columns.')

    psycopg2.extras.execute_batch(cur,
                                  'UPDATE resistome.annotations '
                                  'SET valid = FALSE, valid_reason = %s '
                                  'WHERE annotation_id = %s',
                                  argslist=invalid_annotation_entries,
                                  page_size=10000)

    print('Finished data validation. %i/%i annotations appear to be invalid, with %2.3f passing the validator.'
          % (len(invalid_annotation_entries), count, (1.0 - len(invalid_annotation_entries)/count) * 100.0))
    print('Finished loading data into resistome SQL schema. These entries are marked as VALID = FALSE in the '
          'resistome.annotations table. All other entries appear to be valid.')


if __name__ == '__main__':

    # note: if you are running manually, be sure to update the support DBs first?

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
        cursor.execute('DROP SCHEMA RESISTOME CASCADE')
        with open(os.path.join(constants.INPUT_DIR, 'sql', 'resistome_sql_schema.sql'), 'r') as f:
            sql_schema = ''.join(f.readlines())
            # make schema species specific
            cursor.execute(sql_schema)
        main(cur=cursor)

    connect.commit()
    connect.close()