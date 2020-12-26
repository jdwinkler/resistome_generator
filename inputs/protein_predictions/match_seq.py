import psycopg2
import psycopg2.extras
from resistome import constants
import Bio.SeqIO


if __name__ == '__main__':

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    output_mapping = []

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:

        subfolders = ['MG1655', 'REL606', 'W']
        import glob
        import os

        for strain in subfolders:

            sequence_to_biocyc = dict()

            for filepath in glob.glob(os.path.join(os.getcwd(), strain, strain + '.*')):

                if '.snap2' in filepath:
                    continue

                seqrec = Bio.SeqIO.read(filepath, format='fasta')
                biocyc_accession = seqrec.id.split('|')[-1]

                sequence_to_biocyc[str(seqrec.seq)] = biocyc_accession

            print('Running query for %s' % strain)

            cursor.execute('select genes.accession, aa_seq from genes '
                           'inner join aa_sequences on aa_sequences.gene_id = genes.gene_id '
                           'inner join strain on strain.strain_id = genes.strain_id '
                           'WHERE strain.strain = %s and aa_seq = ANY(%s)',
                           (strain.lower(), list(sequence_to_biocyc.keys())))

            hits = set()
            for record in cursor:

                output_mapping.append((strain,
                                       sequence_to_biocyc[record['aa_seq']],
                                       record['accession']))

                hits.add(record['aa_seq'])

            # missing sequences-need to replace with re-run SNAP2 data for all strains.
            for seq, biocyc in sequence_to_biocyc.items():
                if seq not in hits:
                    print(biocyc, seq)

        with open('biocyc_to_accession_map.txt', 'w') as f:
            for line in output_mapping:
                f.write('\t'.join(line) + '\n')
