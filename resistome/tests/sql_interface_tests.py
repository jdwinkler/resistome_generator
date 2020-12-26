import unittest
from resistome import constants
import psycopg2
import psycopg2.extras
from resistome.sql import sql_interface
import os
from collections import defaultdict


class TestInformationalMethods(unittest.TestCase):

    def test_get_genetic_code(self):

        code, aa_to_codon = sql_interface.get_genetic_code()

        aa_list = constants.AMINO_ACIDS

        # 1 for stop codon

        self.assertTrue(len(set(aa_list) - set(aa_to_codon.keys())) == 1)
        self.assertTrue(len(code.keys()) == 64)

    def test_get_gene_identifiers(self):

        genes = len(sql_interface.get_gene_identifiers())

        self.assertEqual(genes, 4454)

    def test_get_gene_coding_sequence(self):

        test_sequence = 'MVKKSEFERGDIVLVGFDPASGHEQQGAGRPALVLSVQAFNQLGMTLVAPITQGGNFARYAG' \
                        'FSVPLHCEEGDVHGVVLVNQVRMMDLHARLAKRIGLAADEVVEEALLRLQAVVE'

        test_sequence_accession = 'B4225'

        retrieved_sequence = sql_interface.get_polypeptide_sequence(test_sequence_accession)

        self.assertEqual(retrieved_sequence, test_sequence)

    def test_filter_genes_by_feature(self):

        go_tag = 'GO:0016987'

        all_genes = sql_interface.get_gene_identifiers()

        gene_to_go_tuples = sql_interface.get_gene_label_relationships('go')

        genes_with_tag = set([x[0] for x in filter(lambda x: go_tag in x[1] and x[0] in all_genes, gene_to_go_tuples)])

        auto_filtered_genes = sql_interface.filter_genes_by_feature(all_genes,
                                                                    [go_tag],
                                                                    type_of_feature='go')

        self.assertEqual(auto_filtered_genes, genes_with_tag)

    def test_get_gene_interactions(self):

        manual_regulon_dict = defaultdict(list)
        manual_regulon_dict_no_duplicates = defaultdict(set)

        all_genes = sql_interface.get_gene_identifiers()

        with open(os.path.join(constants.INPUT_DIR, 'interactions', 'regulon_db_standardized.txt'), 'r') as f:

            for line in f:
                tokens = line.strip().split('\t')
                regulator = tokens[0]
                regulated = tokens[1]
                direction = tokens[2]

                if regulator not in all_genes:
                    continue

                manual_regulon_dict[regulator].append((regulated, direction))
                manual_regulon_dict_no_duplicates[regulator].add((regulated, direction))

        generated_regulon_dict = sql_interface.get_gene_regulons(all_genes)

        for gene in generated_regulon_dict:

            self.assertTrue(gene in manual_regulon_dict)
            self.assertTrue(gene in manual_regulon_dict_no_duplicates)

            duplicated_element = None
            seen_once = set()

            for element in generated_regulon_dict[gene]:

                if element in seen_once:
                    duplicated_element = element
                seen_once.add(element)

            self.assertEqual(len(generated_regulon_dict[gene]),
                             len(manual_regulon_dict_no_duplicates[gene]),
                             msg='Gene: %s, %s' % (gene, str(duplicated_element)))

            self.assertEqual(set(generated_regulon_dict[gene]),
                             set(manual_regulon_dict_no_duplicates[gene]))

    def test_get_doe_types(self):

        assembled_DOE_types = set()

        with open(os.path.join(constants.INPUT_DIR, 'settings', 'Term Usage.txt'), 'r') as f:
            for line in f:
                tokens = line.strip().split('\t')
                if tokens[0] == 'DOE':
                    assembled_DOE_types.add(tokens[1])

        doe_types = sql_interface.get_paper_doe_types()

        self.assertEqual(set(doe_types), assembled_DOE_types)

    def test_get_phenotype_classes(self):

        assembled_phenotype_classes = set()

        with open(os.path.join(constants.INPUT_DIR, 'settings', 'Term Usage.txt'), 'r') as f:
            for line in f:
                tokens = line.strip().split('\t')
                if tokens[0] == 'Tags':
                    assembled_phenotype_classes.add(tokens[1])

        phenotype_classes = sql_interface.get_phenotype_classes()

        self.assertEqual(set(phenotype_classes), assembled_phenotype_classes)


class TestDataExtractionMethods(unittest.TestCase):

    def test_get_mutant_genotypes(self):

        database_name = constants.DB_NAME
        user_name = constants.DB_USERNAME
        password = constants.DB_PASSWORD

        connection = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'"
                                      % (database_name, user_name, password))
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        cursor.execute('select mutant_id from resistome.mutants')
        assembled_mutant_ids = [x['mutant_id'] for x in cursor.fetchall()]

        data = sql_interface.get_mutant_genotypes()
        mutant_ids = [x['mutant_id'] for x in data]

        self.assertEqual(sorted(assembled_mutant_ids), sorted(mutant_ids))

        assembled_genotypes = dict()

        for mutant_id in assembled_mutant_ids:

            cursor.execute('select * from resistome.mutations where mutant_id = %s', (mutant_id,))
            mutations = cursor.fetchall()

            gene_ids = [x['gene_id'] for x in mutations]
            genes = [x['name'] for x in mutations]

            species_agg = []
            mutation_types = []
            annotations = []

            for gene_id, gene in zip(gene_ids, genes):
                cursor.execute('select * from resistome.annotations where gene_id = %s', (gene_id,))
                gene_data = cursor.fetchall()

                species_agg.extend([gene] * len(gene_data))
                mutation_types.extend([x['mutation_type'] for x in gene_data])
                annotations.extend([x['annotation'] for x in gene_data])

            cursor.execute('select * from resistome.phenotypes where mutant_id = %s', (mutant_id,))
            phenotype_data = cursor.fetchall()

            phenotypes = [x['phenotype'] for x in phenotype_data]
            pclasses = [x['phenotype_class'] for x in phenotype_data]
            ptypes = [x['phenotype_type'] for x in phenotype_data]

            assembled_genotypes[mutant_id] = {'mutations': mutation_types,
                                              'annotations': annotations,
                                              'species_names': species_agg,
                                              'phenotypes': phenotypes,
                                              'pclasses': pclasses,
                                              'types': ptypes}
        for genotype in data:

            mutant_id = genotype['mutant_id']

            self.assertEqual(sorted(genotype['mutations']),
                             sorted(assembled_genotypes[mutant_id]['mutations']))

            self.assertEqual(sorted(genotype['annotations']),
                             sorted(assembled_genotypes[mutant_id]['annotations']))

            self.assertEqual(sorted(genotype['paired_species_names']),
                             sorted(assembled_genotypes[mutant_id]['species_names']))

            self.assertEqual(sorted(genotype['phenotypes']),
                             sorted(assembled_genotypes[mutant_id]['phenotypes']))

            self.assertEqual(sorted(genotype['pclasses']),
                             sorted(assembled_genotypes[mutant_id]['pclasses']))

            self.assertEqual(sorted(genotype['types']),
                             sorted(assembled_genotypes[mutant_id]['types']))

    def test_get_gene_mutation_tuples(self):

        database_name = constants.DB_NAME
        user_name = constants.DB_USERNAME
        password = constants.DB_PASSWORD

        connection = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'"
                                      % (database_name, user_name, password))
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        cursor.execute('select annotation_id from resistome.annotations')
        assembled_annotation_ids = [x['annotation_id'] for x in cursor.fetchall()]
        record_count = len(assembled_annotation_ids)

        data = sql_interface.get_gene_mutation_tuples()

        self.assertEqual(len(data), record_count)

        assembled_annotation_dict = defaultdict(list)
        annotation_dict = defaultdict(list)

        for annotation_id in assembled_annotation_ids:
            cursor.execute('select gene_id, mutation_type, annotation from resistome.annotations where annotation_id = %s', (annotation_id,))
            mutation_tuple = cursor.fetchall()
            for result in mutation_tuple:
                assembled_annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))

        for result in data:
            annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))
        self.assertEqual(sorted(annotation_dict.keys()), sorted(assembled_annotation_dict.keys()))
        for gene_id in annotation_dict:
            self.assertEqual(sorted(annotation_dict[gene_id]), sorted(assembled_annotation_dict[gene_id]))

        mutation_of_interest = ['aa_snps']
        filtered_assembled_annotation_dict = defaultdict(list)

        for key in assembled_annotation_dict:
            for data_tuple in assembled_annotation_dict[key]:
                if data_tuple[0] == 'aa_snps':
                    filtered_assembled_annotation_dict[key].append(data_tuple)

        data = sql_interface.get_gene_mutation_tuples(filter_for_mutation_types=mutation_of_interest)
        annotation_dict = defaultdict(list)

        for result in data:
            annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))
        self.assertEqual(sorted(annotation_dict.keys()), sorted(filtered_assembled_annotation_dict.keys()))
        for gene_id in annotation_dict:
            self.assertEqual(sorted(annotation_dict[gene_id]), sorted(filtered_assembled_annotation_dict[gene_id]))

        phenotype_class = 'solvents_biofuels'

        annotation_ids_to_parse = set()

        for annotation_id in assembled_annotation_ids:

            cursor.execute('select gene_id from resistome.annotations where annotation_id = %s', (annotation_id,))
            gid_results = cursor.fetchall()
            for result in gid_results:
                gene_id = result['gene_id']
                cursor.execute('select mutant_id from resistome.mutations where gene_id = %s', (gene_id, ))

                mid_results = cursor.fetchall()
                for mid_result in mid_results:
                    cursor.execute('select array_agg(phenotype_class) as pclasses from resistome.phenotypes where mutant_id = %s group by mutant_id', (mid_result['mutant_id'],))
                    phenotype_results = cursor.fetchone()

                    if phenotype_class in phenotype_results['pclasses']:
                        annotation_ids_to_parse.add(annotation_id)

        assembled_annotation_dict = defaultdict(list)
        for annotation_id in annotation_ids_to_parse:
            cursor.execute('select gene_id, mutation_type, annotation from resistome.annotations where annotation_id = %s', (annotation_id,))
            mutation_tuple = cursor.fetchall()
            for result in mutation_tuple:
                assembled_annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))

        data = sql_interface.get_gene_mutation_tuples(filter_for_phenotype_class=set([phenotype_class]))
        annotation_dict = defaultdict(list)

        for result in data:
            annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))
        self.assertEqual(sorted(annotation_dict.keys()), sorted(assembled_annotation_dict.keys()))
        for gene_id in annotation_dict:
            self.assertEqual(sorted(annotation_dict[gene_id]), sorted(assembled_annotation_dict[gene_id]))

        assembled_annotation_dict = defaultdict(list)
        for annotation_id in annotation_ids_to_parse:
            cursor.execute('select gene_id, mutation_type, annotation from resistome.annotations where annotation_id = %s', (annotation_id,))
            mutation_tuple = cursor.fetchall()
            for result in mutation_tuple:
                if result['mutation_type'] == 'aa_snps':
                    assembled_annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))

        data = sql_interface.get_gene_mutation_tuples(filter_for_phenotype_class=set([phenotype_class]),
                                                      filter_for_mutation_types=set(['aa_snps']))
        annotation_dict = defaultdict(list)

        for result in data:
            annotation_dict[result['gene_id']].append((result['mutation_type'], result['annotation']))
        self.assertEqual(sorted(annotation_dict.keys()), sorted(assembled_annotation_dict.keys()))
        for gene_id in annotation_dict:
            self.assertEqual(sorted(annotation_dict[gene_id]), sorted(assembled_annotation_dict[gene_id]))

    def test_get_mutation_location_data(self):

        database_name = constants.DB_NAME
        user_name = constants.DB_USERNAME
        password = constants.DB_PASSWORD

        connection = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'"
                                      % (database_name, user_name, password))
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        location_data = dict()

        for species in constants.SPECIES_LIST:

            temp_data = sql_interface.get_gene_location_data(query_species=species)

            for gene in temp_data:

                location_data[gene] = ('chr1',
                                     temp_data[gene]['start'],
                                     temp_data[gene]['stop'])

        mutation_data = sql_interface.get_mutation_location_data(location_data, debug_flag=True)

        cursor.execute('select name, gene_id from resistome.mutations')
        gene_data = cursor.fetchall()

        assembled_gene_dict = defaultdict(list)

        for result in gene_data:

            gene_name = result['name']
            gene_id = result['gene_id']

            cursor.execute('select annotation, mutation_type from resistome.annotations where gene_id = %s', (gene_id,))
            m_a_data = cursor.fetchall()

            for ma_record in m_a_data:
                assembled_gene_dict[gene_name].append((ma_record['annotation'],
                                             ma_record['mutation_type']))

        gene_dict = defaultdict(list)

        for result in mutation_data:

            for annotation, mutation_type in zip(result['annotations'], result['mutation_types']):
                gene_dict[result['name']].append((annotation, mutation_type))

        self.assertTrue(len(set(gene_dict.keys()) - set(assembled_gene_dict.keys())) == 0)

        for gene in gene_dict:

            self.assertEqual(sorted(gene_dict[gene]), sorted(assembled_gene_dict[gene]))

    def test_mutation_by_doe_type(self):

        database_name = constants.DB_NAME
        user_name = constants.DB_USERNAME
        password = constants.DB_PASSWORD

        connection = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'"
                                      % (database_name, user_name, password))
        cursor = connection.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        classes_to_process = ['none', 'modeled', 'statistical', 'completionist', 'parallel', 'serial', 'random']

        for design_class in classes_to_process:

            assembled_mutation_type_counts = []
            assembled_gene_counts = []

            cursor.execute('select paper_id from resistome.papers where methods <@ %s', ([design_class],))
            paper_ids = [x['paper_id'] for x in cursor]

            cursor.execute('select mutant_id from resistome.mutants where paper_id = ANY(%s)', (paper_ids,))
            mutant_ids = [x['mutant_id'] for x in cursor]

            cursor.execute('select '
                           'array_agg(distinct resistome.mutations.gene_id) as gene_ids,'
                           'array_agg(distinct resistome.mutations.name) as species_names, '
                           'array_agg(resistome.annotations.mutation_type) as mutation_types '
                           'from resistome.mutations '
                           'inner join resistome.annotations on (resistome.annotations.gene_id = resistome.mutations.gene_id) '
                           'where mutant_id = ANY(%s) '
                           'group by resistome.mutations.mutant_id',
                           (mutant_ids,))

            for result in cursor:
                assembled_mutation_type_counts.append(len(result['mutation_types']))
                assembled_gene_counts.append(len(result['species_names']))

            mutation_type_counts, gene_counts = sql_interface.get_mutation_counts_by_doe_type(types_to_include=[design_class])

            self.assertEqual(len(mutation_type_counts), len(assembled_mutation_type_counts))
            self.assertEqual(len(gene_counts), len(assembled_gene_counts))

            self.assertEqual(sorted(mutation_type_counts), sorted(assembled_mutation_type_counts))
            self.assertEqual(sorted(gene_counts), sorted(assembled_gene_counts))




