from resistome.sql import biocyc_data_parser, resistome_builder, ncbi_data_parser, public_builder
from resistome.utils import name_helper
from resistome import constants
import os
import psycopg2
import psycopg2.extras

"""

This script automates construction of the support and Resistome databases. See readme.md for more details.

"""

if __name__ == '__main__':

    print('Using the following credentials to create the Resistome database:')
    print('Credential file path: %s' % os.path.join(constants.INPUT_DIR, 'db_credentials'))
    print('Building Biocyc/NCBI databases')

    root_dir = constants.INPUT_DIR
    ECOLI_MG1655_DIR = os.path.join(root_dir, 'biocyc', 'mg1655')
    ECOLI_B_DIR = os.path.join(root_dir, 'biocyc', 'rel606')

    ECOLI_BW25113_DIR = os.path.join(root_dir, 'biocyc', 'bw25113')
    ECOLI_W3110_DIR = os.path.join(root_dir, 'biocyc', 'w3110')
    ECOLI_MDS42_DIR = os.path.join(root_dir, 'biocyc', 'mds42')
    ECOLI_BL21_DIR = os.path.join(root_dir, 'biocyc', 'bl21')
    ECOLI_BL21_DE3_DIR = os.path.join(root_dir, 'biocyc', 'bl21_de3_')

    source_dirs = [ECOLI_W3110_DIR, ECOLI_BW25113_DIR, ECOLI_MDS42_DIR, ECOLI_BL21_DIR, ECOLI_BL21_DE3_DIR,
                   ECOLI_MG1655_DIR, ECOLI_B_DIR]
    # source_dirs = [ECOLI_W3110_DIR]
    ncbi_strain_tuples = [(x, os.path.split(x)[1]) for x in source_dirs]

    # internally we use W for this strain alias, unlike the others.
    ncbi_strain_tuples.append((os.path.join(root_dir, 'biocyc', 'ecoliw'), 'w'))

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:

        cursor.execute('DROP SCHEMA IF EXISTS PUBLIC CASCADE')
        cursor.execute('DROP SCHEMA IF EXISTS RESISTOME CASCADE')

        schemas_to_create = ['resistome_biocyc_public.sql',
                             'resistome_biocyc_schema.sql',
                             'resistome_sql_schema.sql']

        for schema_file in schemas_to_create:
            with open(os.path.join(constants.INPUT_DIR, 'sql', schema_file)) as f:
                sql_schema = ''.join(f.readlines())
                cursor.execute(sql_schema)

        try:
            public_builder.main(cursor)
            for source_dir, strain in ncbi_strain_tuples:
                print('Processing NCBI annotations for: %s' % strain)
                ncbi_data_parser.main(cur=cursor, source_data=source_dir, species='Escherichia coli', strain=strain)
        except:
            print('Exception occurred when generating biocyc supporting databases.')
            print('Possible causes include missing data (biocyc data is distributed separately from the rest of the database '
                  'for copyright reasons), the database credentials are wrong, the host database has not been set up, or '
                  'an unknown cause.')
            raise

        print('Succeeded building Biocyc/NCBI database!')

        print('\nGenerating name mapping tables to B-numbers...')

        try:
            name_helper.main(cursor)
        except:
            print('Failed when generating name => accession files!')
            raise

        print('Name mapping generated successfully!')

        try:
            print('\nGenerating Sauer Metabolite tables...')
            public_builder.build_metabolite_table(cursor)
            print('Succeeded in building Sauer Metabolite tables!')
        except:
            print('Failed when generating metabolite tables!')
            raise

        try:
            print('\nBuilding Resistome tables')
            resistome_builder.main(cursor, add_resistome_go_metabolite_tables=True)
        except:
            print('Exception occurred when building Resistome database')
            print('This is most likely caused by the target Postgres database not being set up, or new data failing QC.')
            raise

    connect.commit()
    connect.close()
