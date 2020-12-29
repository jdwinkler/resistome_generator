from resistome.utils import database_utils
import os
import psycopg2
import psycopg2.extras
from resistome import constants
from typing import List, Tuple
from collections import defaultdict


def validate_mutation_data(cursor) -> Tuple[List[Tuple[int, str]], int]:

    """

    Checks the validity for mutations the are specified using verifable information, such as gene locations, bases,
    or residues. Returns a List of tuples (reason why, invalid annotation_id pk) meant to use in an update query.

    You can see below to determine how the valid/invalid calls are made.

    :param cursor:
    :return:
    """

    cursor.execute('select resistome.annotations.*, '
                   'resistome.mutations.name, '
                   'resistome.mutations.strain '
                   'FROM resistome.mutations '
                   'INNER JOIN resistome.annotations on resistome.annotations.gene_id = resistome.mutations.gene_id')

    mutation_data = cursor.fetchall()

    strains = set()

    for record in mutation_data:
        strains.add(record['strain'])

    strains = list(strains)

    aa_sequences = defaultdict(dict)
    genome_sequences = dict()
    gene_locations = defaultdict(dict)

    for strain in strains:

        cursor.execute('select sequence, strain.strain from genome '
                       'INNER JOIN strain on genome.strain_id = strain.strain_id '
                       'WHERE strain.strain = lower(%s)', (strain, ))
        genome_sequences[strain] = cursor.fetchone()['sequence']

        cursor.execute('select genes.accession, aa_sequences.aa_seq '
                       'FROM genes '
                       'INNER JOIN aa_sequences on aa_sequences.gene_id = genes.gene_id '
                       'INNER JOIN strain on genes.strain_id = strain.strain_id '
                       'WHERE strain.strain = lower(%s)', (strain, ))

        for record in cursor:
            aa_sequences[strain][record['accession']] = record['aa_seq']

        cursor.execute('select genes.accession, gene_locations.start, gene_locations.stop '
                       'FROM genes '
                       'INNER JOIN gene_locations on gene_locations.gene_id = genes.gene_id '
                       'INNER JOIN strain on genes.strain_id = strain.strain_id '
                       'WHERE strain.strain = lower(%s)', (strain,))

        for record in cursor:
            gene_locations[strain][record['accession']] = (record['start'], record['stop'])

    invalid_annotations = []

    missing_gene_counter = defaultdict(int)

    for record in mutation_data:

        annotation_id = record['annotation_id']
        strain = record['strain']
        gene_name = record['name']
        mutation_type = record['mutation_type']
        annotation_dict = record['annotation']

        if gene_name not in gene_locations[strain]:
            # can't find the correct gene
            missing_gene_counter[(strain, gene_name)] += 1
            invalid_annotations.append((annotation_id,  'MISSING GENE'))
            continue

        if mutation_type == 'aa_snps':
            for (location, wt_aa, alt_aa) in annotation_dict[mutation_type]:
                aa_seq = aa_sequences[strain].get(gene_name, None)
                if aa_seq is None:
                    # no sequence
                    invalid_annotations.append((annotation_id, 'MISSING AA SEQ'))
                else:
                    if len(aa_seq) < location:
                        # aa seq is shorter than the residue position specified
                        invalid_annotations.append((annotation_id, 'AA POSITION NOT IN PROTEIN'))
                    elif len(aa_seq) == location:
                        if wt_aa != '*':
                            # codon specified is last but data claims that this is not the stop codon
                            invalid_annotations.append((annotation_id, 'AA POSITION LAST BUT NOT STOP'))
                    else:
                        actual_wt_aa = aa_seq[location]
                        if actual_wt_aa != wt_aa:
                            # claimed/actual WT residue are different
                            invalid_annotations.append((annotation_id, 'WT/REF AA MISMATCH'))
        if mutation_type == 'nuc_snps':
            for (location, wt_base, alt_base) in annotation_dict[mutation_type]:
                if location is None:
                    # no location given
                    invalid_annotations.append((annotation_id, 'NO SNP LOCATION'))
                else:
                    if location < 0 or location > len(genome_sequences[strain]):
                        # location does not fall in genome
                        invalid_annotations.append((annotation_id, 'SNP POSITION NOT IN GENOME'))
                    else:
                        actual_wt_base = genome_sequences[strain][location]
                        if actual_wt_base != wt_base:
                            # base mismatch between claimed/actual WT base
                            invalid_annotations.append((annotation_id, 'WT/REF BASE MISMATCH'))
        if mutation_type == 'indel':
            for entry in annotation_dict[mutation_type]:
                if entry[-1] == 'gene_unknown':
                    # can't find gene
                    invalid_annotations.append((annotation_id, 'INDEL IN UNKNOWN GENE'))
                elif entry[0] is None:
                    # no location provided
                    invalid_annotations.append((annotation_id, 'INDEL IN UNKNOWN LOCATION'))
                elif entry[0] < 0 or entry[0] > len(genome_sequences[strain]):
                    # indel not in genome
                    invalid_annotations.append((annotation_id, 'INDEL NOT IN GENOME'))
        if mutation_type == 'duplication':
            for entry in annotation_dict[mutation_type]:
                # see above for indel explanation, basically the same here
                if entry[0] < 0 or entry[0] > len(genome_sequences[strain]):
                    invalid_annotations.append((annotation_id, 'DUPLICATION LEFTPOS NOT IN GENOME'))
                if entry[1] < 0 or entry[0] > len(genome_sequences[strain]):
                    invalid_annotations.append((annotation_id, 'DUPLICATION RIGHTPOS NOT IN GENOME'))
        if 'large_' in mutation_type:
            for gene in annotation_dict[mutation_type]:
                if isinstance(gene, int):
                    # skip, last entry of large_amplification is supposed to fold amp
                    break
                if gene not in gene_locations[strain]:
                    # missing gene specified in large mutation definition
                    invalid_annotations.append((annotation_id, '%s MISSING GENE' % mutation_type.upper()))
                    break
            if mutation_type == 'large_amplification' and not isinstance(annotation_dict[mutation_type][-1], int):
                # missing fold X entry for large amplification
                invalid_annotations.append((annotation_id, 'LARGE AMP MISSING X-FOLD ENTRY'))

    invalid_annotations = [(y, x) for x, y in invalid_annotations]

    for x in sorted(missing_gene_counter, key=lambda x: missing_gene_counter[x], reverse=True):
        print(x, missing_gene_counter[x])

    return invalid_annotations, len(mutation_data)


if __name__ == '__main__':

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    cur = connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    errors, count = validate_mutation_data(cur)

    for x in errors:
        print(x)

    print('Number of errors total: %i/%i' % (len(errors), count))
