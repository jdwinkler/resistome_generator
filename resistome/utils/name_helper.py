import psycopg2
import psycopg2.extras
from collections import defaultdict
from typing import Dict, Tuple, Set
import os

from resistome import constants


def add_is_hyphens(gene_name: str):

    """

    Generates hyphen free and hyphenated versions of IS element names to improve gene capture.

    :param gene_name:
    :return:
    """

    output_names = set()
    if str.isnumeric(gene_name[-1]) and 'INS' in gene_name.upper():

        if '-' in gene_name:
            output_names.add(gene_name.replace('-', ''))
        else:
            new_name = []
            alpha_to_number = False
            for c in gene_name:
                if str.isnumeric(c) and not alpha_to_number:
                    new_name.append('-')
                    alpha_to_number = True
                new_name.append(c)
            new_name = ''.join(new_name)
            output_names.add(new_name)

    return output_names


def get_mg1655_accessions(cursor, targets, source_to_accession, already_blacklisted):

    """

    Tries to map each gene accession to a unique MG1655 bnumber. Mappings are only made if a single unique mapping
    exists, or can be inferred from the accession (when using BW25113 accessions) or ECK unique identifier.

    If the synonyms for a given accession already include a B-number, that is used instead of querying the database.

    :param cursor:
    :param targets:
    :param source_to_accession:
    :param already_blacklisted:
    :return:
    """

    accessions_to_synonym = defaultdict(set)
    for x, y in source_to_accession.items():
        accessions_to_synonym[y].add(x)

    collated_mapping = defaultdict(set)
    for accession, synonyms in accessions_to_synonym.items():

        synonyms = sorted([x.upper() for x in synonyms])

        bnumber = None
        for syn in synonyms:
            # get the bnumber
            if len(syn) == 5 and syn[0].upper() == 'B':
                bnumber = syn
                break

        if bnumber is not None:
            # if one exists, just map all the synonyms to it then go onto the next gene
            for syn in synonyms:
                collated_mapping[syn].add(bnumber)
            continue

        # select all the E. coli MG1655 genes that have overlapping synoyms, the same gene name, or the same accession
        cursor.execute('select name, synonyms, accession, start from genes '
                       'inner join gene_locations on gene_locations.gene_id = genes.gene_id '
                       'WHERE strain_id = ANY(select strain_id from strain where strain.strain = %s) AND '
                       '(genes.synonyms @> %s or genes.name = ANY(%s) or genes.accession = ANY(%s))',
                       ('mg1655', synonyms, synonyms, synonyms))

        result = cursor.fetchall()

        if len(result) > 1:
            # Non-unique; try to match if possible.
            eck = list(filter(lambda x: 'ECK' in x, synonyms))
            eck_matches = []

            for rec in result:
                # print('non-unique mapping', rec)
                if len(eck) == 1 and eck[0] in rec['synonyms']:
                    eck_matches.append(rec)
            if len(eck_matches) == 1:
                # ECK attempt
                print('matched based on ECK number: %s' % accession)
                for x in synonyms:
                    collated_mapping[x].add(eck_matches[0]['accession'])
            elif len(eck_matches) > 0 and 'BW25113' in accession:
                # BW25113 accessions follow this pattern (although it is sometimes wrong)
                guessed_bnumber = 'B' + accession.split('_')[-1]
                for rec in eck_matches:
                    if rec['accession'] == guessed_bnumber:
                        print('matched based on inferred B-number: %s for %s' % (guessed_bnumber, accession))
                        for x in synonyms:
                            collated_mapping[x].add(rec['accession'])
            elif len(synonyms) == 1:
                # try an exact name match iff there is a single name.
                for rec in result:
                    if rec['name'] == synonyms[0]:
                        print('matched based on exact name match')
                        collated_mapping[rec['name']].add(rec['accession'])
            elif 'TSAA' in synonyms:
                print('manual match for TSAA')
                collated_mapping['TSAA'].add('B0195')
            else:
                print('skipped rec, no match inference strategy: %s, %s' % (accession, synonyms))
                # for x in result:
                #     print(accession, synonyms, x)
        elif len(result) > 0:
            for x in synonyms:
                collated_mapping[x].add(result[0]['accession'])

    output = dict()
    for x, y in collated_mapping.items():

        if len(y) > 1:
            print(x, y, 'mg1655')
            raise AssertionError('Did not obtain unique mapping for MG1655 synonyms!')
        else:
            output[x] = y.pop()

    # check to make sure the cross-species mappings are unique
    accession_to_bnumber_check = defaultdict(set)
    for x, y in output.items():
        if x in source_to_accession:
            accession_to_bnumber_check[source_to_accession[x]].add(y)

    fail = False
    for x, y in accession_to_bnumber_check.items():
        if len(y) > 1:
            print(x, y, 'accession mapping')
            fail = True

    if fail:
        raise AssertionError('Failed to map accessions uniquely!')

    return output


def extract_name_mapping(strain: str, cursor):
    """

    Extracts name mappings to help disambiguate input gene names. Pulls data from corresponding species tables
    after running biocyc/NCBI data loads.

    Tries to add synonyms (if they do not overlap with official names) and common variants of insertion sequence
    names.

    :param strain:
    :param connect:
    :return:
    """

    cursor.execute('select name, synonyms, accession, start from genes '
                   'inner join gene_locations on gene_locations.gene_id = genes.gene_id '
                   'INNER JOIN strain on strain.strain_id = genes.strain_id '
                   'WHERE strain.strain = %s', (strain,))

    records = cursor.fetchall()
    mapping_table = defaultdict(set)

    eck_jw_to_original_accession = dict()

    blacklist_synonyms = blacklist_nonspecific_synonyms(records)
    accession_to_start = dict()

    if strain == 'rel606':
        # manually map this...
        mapping_table['YLCI'].add('ECB_00510')
        accession_to_start['YLCI'] = 550723

    for record in records:

        # print(record)
        mapping_table[record['name']].add(record['accession'])
        accession_to_start[record['accession']] = record['start']

        if strain != 'mg1655':
            found_bnumber = False

            # try to make this deterministic
            for synonym in sorted(record['synonyms']):

                if synonym.upper() in blacklist_synonyms:
                    continue

                mapping_table[synonym].add(record['accession'])
                try:
                    if synonym[0] == 'B' and len(synonym) == 5:
                        # TODO make this a function
                        found_bnumber = True
                except IndexError:
                    print(synonym, strain, record)
                    raise

            if not found_bnumber:
                # if we cannot extract the bumber from the synonyms, try using the ECK or JW accessions instead.
                eck_jw_number = filter(lambda x: 'JW' in x or 'ECK' in x, record['synonyms'])
                for x in eck_jw_number:
                    eck_jw_to_original_accession[x] = record['accession']
        else:
            for synonym in sorted(record['synonyms']):
                if synonym.upper() in blacklist_synonyms:
                    # ignore any non-specific synonyms.
                    continue
                mapping_table[synonym].add(record['accession'])

            if 'TSAA' in record['synonyms'] and strain == 'mg1655' and record['name'] == 'TRMO':
                # manual matching since this is often mutated?
                mapping_table['TSAA'].add(record['accession'])

    if len(eck_jw_to_original_accession.keys()) > 0:
        # cross reference using these alternative unique IDs (ECK is ecocyc, JW is the Keio collection)
        cursor.execute('select name, synonyms, accession, start from genes '
                       'inner join gene_locations on gene_locations.gene_id = genes.gene_id '
                       'WHERE strain_id = ANY(select strain_id from strain where strain.strain = %s) AND '
                       'genes.synonyms @> %s', ('mg1655', list(eck_jw_to_original_accession)))

        for record in cursor:
            for syn in sorted(record['synonyms']):
                if 'ECK' in syn or 'JW' in syn:
                    # mapping_table[record['name']].add(eck_jw_to_original_accession[syn])
                    mapping_table[record['accession']].add(eck_jw_to_original_accession[syn])

    accessions_to_update = dict()
    for x, y in mapping_table.items():
        if len(y) > 1:
            if strain == 'w':
                # ignore plasmid annotations if they conflict with the chromosomal ones.
                removed_plasmid_accession = list(filter(lambda z: 'ECW_P' not in z, y))
                if len(removed_plasmid_accession) != 1:
                    raise AssertionError(y)
                accessions_to_update[x] = removed_plasmid_accession[0]
            elif 'INS' not in x:
                # map the accessions to themselves alone, remove all synonyms
                for loci in y:
                    accessions_to_update[loci] = loci
                    accessions_to_update[x] = None
            elif strain == 'bw25113':
                # rename IS elements to something more logical
                ordered_loci = sorted(y, key=lambda x: int(x.split('_')[1]))
                accessions_to_update[x] = ordered_loci[0]

                if not str.isnumeric(x[-1]):
                    accessions_to_update[x + '1'] = ordered_loci[0]

                k = 2
                while k - 1 < len(ordered_loci):
                    accession = ordered_loci[k - 1]
                    new_name = x[0:4] + str(k)
                    accessions_to_update[new_name] = accession
                    k += 1

    nonspecific_genes = set()
    nonspecific_genes.update(blacklist_synonyms)
    for x, y in accessions_to_update.items():
        if y is None:
            # remove anything masked off in the code above
            del mapping_table[x]
            nonspecific_genes.add(x)
        else:
            mapping_table[x] = {y}

    output_table = dict()
    for x, y in mapping_table.items():
        if len(y) != 1:
            # check to make sure each gene is only mapped to a single accession
            print(x, y, '*****')
            raise AssertionError('Failed to handle case of gene => {id 1, id 2} multiple mapping.')
        else:
            output_table[x] = y.pop()
            if 'INS' in x:
                # add hyphens to capture common name variants
                for is_insert in add_is_hyphens(x):
                    output_table[is_insert] = output_table[x]

    mapping_to_mg1655 = get_mg1655_accessions(cursor,
                                              source_to_accession=output_table,
                                              targets=list(output_table.keys()),
                                              already_blacklisted=blacklist_synonyms)

    return output_table, mapping_to_mg1655, nonspecific_genes, accession_to_start


def blacklist_nonspecific_synonyms(records):

    """

    Blacklist any common synonyms that map to multiple genes within an assembly.

    :param records:
    :return:
    """

    synonym_to_accession_tester = defaultdict(set)
    blacklist_synonyms = set()
    for record in records:
        for synonym in record['synonyms']:
            synonym_to_accession_tester[synonym.upper()].add(record['accession'])
            synonym_to_accession_tester[synonym.upper().strip()].add(record['accession'])
        synonym_to_accession_tester[record['name'].upper()].add(record['accession'])
        synonym_to_accession_tester[record['name'].upper().strip()].add(record['accession'])
    for x, y in synonym_to_accession_tester.items():
        # any gene synonym that maps to multiple accession in the genome.
        if len(y) > 1:
            blacklist_synonyms.add(x.upper())

    return blacklist_synonyms


def update_with_mg1655_names(mapping_table, to_mg1655, blacklist, cursor):

    """

    Backpropagate mg1655 names onto the target species-since MG1655 is the best annotation for E. coli, this
    increases the number of gene-accession relationships that we can easily capture.

    :param mapping_table:
    :param to_mg1655:
    :param blacklist:
    :param cursor:
    :return:
    """

    all_bnumbers = list(set(to_mg1655.values()))

    output_dict = dict(mapping_table)
    to_mg1655 = dict(to_mg1655)

    cursor.execute('select accession, name, synonyms from genes '
                   'where accession = any(%s)', (all_bnumbers,))

    bnumber_to_details = dict()
    for record in cursor:
        bnumber_to_details[record['accession']] = (record['name'], record['synonyms'])

    for name, accession in mapping_table.items():

        if name not in to_mg1655:
            # no bnumber mapping original, gotta bail out.
            continue

        mapped_bnumber = to_mg1655[name]

        if mapped_bnumber not in bnumber_to_details:
            # bnumber is not known (?)
            continue

        mg1655_name, mg1655_synonyms = bnumber_to_details[mapped_bnumber]
        synonyms_to_check = {mg1655_name}
        synonyms_to_check.update(mg1655_synonyms)

        for syn_check in synonyms_to_check:

            if syn_check in blacklist:
                continue
            if syn_check in output_dict:
                continue

            # add new outputs.
            output_dict[syn_check] = accession
            to_mg1655[syn_check] = mapped_bnumber

    return output_dict, to_mg1655


def error_check(output_tuples):
    duplicate_gene_strain_ids = set()
    errors = []

    for x in output_tuples:

        if (x[0], x[2]) in duplicate_gene_strain_ids:
            errors.append(x)

        duplicate_gene_strain_ids.add((x[0], x[2]))

    if len(errors) > 0:
        print('Errors detected!')
        for y in errors:
            print(y)
        raise AssertionError('Detected duplicate entries for some gene/strain combinations. This means the same '
                             'identifier is being mapped to multiple different values, which cannot be correct.')


def main(cursor):
    strains = constants.SPECIES_LIST
    species = 'Escherichia coli'

    if 'mg1655' not in strains:
        raise AssertionError('Expected mg1655 to be present in the constants species list, but it is missing?')

    species_rows = []
    sequence_rows = []
    for strain in strains:

        mapping_table, to_mg1655, blacklist, start_locations = extract_name_mapping(strain, cursor)

        if strain != 'mg1655':
            old_count = len(mapping_table.keys())
            mapping_table, to_mg1655 = update_with_mg1655_names(mapping_table, to_mg1655, blacklist, cursor)
            print('Added %i elements to mapping table' % (len(mapping_table.keys()) - old_count,))

        duplicates = set()
        for gene, accession in sorted(mapping_table.items(), key=lambda x: x[0]):
            if (gene.upper().strip(), strain, accession.upper()) not in duplicates:
                duplicates.add((gene.upper().strip(), strain, accession.upper()))
                row = [gene.upper().strip(), 'Escherichia coli', strain.lower(), accession.upper(), start_locations[accession]]
                species_rows.append(row)

            if strain != 'mg1655' and gene in to_mg1655 and (to_mg1655[gene], strain, accession.upper()) not in duplicates:
                if to_mg1655[gene] not in mapping_table:
                    duplicates.add((to_mg1655[gene], strain, accession.upper()))
                    row = [to_mg1655[gene], 'Escherichia coli', strain.lower(), accession.upper(),
                           start_locations[accession]]
                    species_rows.append(row)

        accession_bnumber_pairs = set()
        duplicates = set()
        for gene, accession in sorted(mapping_table.items(), key=lambda x: x[0]):
            if gene in to_mg1655:
                if (gene.upper().strip(), strain, to_mg1655[gene].upper()) not in duplicates:
                    duplicates.add((gene.upper().strip(), strain, to_mg1655[gene].upper()))
                    row = [gene.upper().strip(), species, strain.lower(), to_mg1655[gene], start_locations[accession]]
                    accession_bnumber_pairs.add((accession, to_mg1655[gene]))
                    sequence_rows.append(row)

        for (accession, bnumber) in sorted(accession_bnumber_pairs, key=lambda x: x[0]):
            # accession to bnum
            if (accession.upper(), strain, bnumber.upper()) not in duplicates and accession not in mapping_table:
                duplicates.add((accession.upper(), strain, bnumber.upper()))
                row = [accession.upper().strip(), species, strain.lower(), bnumber, start_locations[accession]]
                sequence_rows.append(row)

            # self => self mapping
            if strain != 'mg1655':
                if (bnumber.upper(), strain, bnumber.upper()) not in duplicates and bnumber not in mapping_table:
                    duplicates.add((bnumber.upper(), strain, bnumber.upper()))
                    row = [bnumber, species, strain.lower(), bnumber, start_locations[accession]]
                    sequence_rows.append(row)

    error_check(species_rows)
    error_check(sequence_rows)

    species_rows = sorted(species_rows, key=lambda x: (x[0], x[2]))
    sequence_rows = sorted(sequence_rows, key=lambda x: (x[0], x[2]))

    with open(os.path.join(constants.INPUT_DIR, 'standardization', 'Species-Specific Gene Names_DB.txt'), 'w') as f:
        for line in species_rows:
            f.write('\t'.join(map(str, line)) + '\n')

    with open(os.path.join(constants.INPUT_DIR, 'standardization', 'Sequence-Specific Gene Names_DB.txt'), 'w') as f:
        for line in sequence_rows:
            f.write('\t'.join(map(str, line)) + '\n')


if __name__ == '__main__':

    try:
        connect = psycopg2.connect("dbname='%s' user='%s' host='localhost' password='%s'" % (constants.DB_NAME,
                                                                                             constants.DB_USERNAME,
                                                                                             constants.DB_PASSWORD))
    except:
        raise

    with connect.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
        main(cursor)

