from typing import List, Dict
from collections import defaultdict
import resistome.sql.sql_interface as sql_interface


def hamming_distance(str1, str2):
    """

    Calculates hamming distance between two different strings of the same length.

    :param str1:
    :param str2:
    :return:
    """

    distance = 0
    if len(str1) != len(str2):
        raise AssertionError('unequal length in strings')

    # calculate hamming distance
    for base_x, base_y in zip(str1, str2):
        if base_x != base_y:
            distance += 1

    return distance


def infer_minimal_snps_from_aa(nt_sequence: str,
                               gene_start: int,
                               wt_residue: str, position: int, alt_aa: str,
                               aa_to_codons: Dict[str, List[str]],
                               codon_to_aa: Dict[str, str]):
    """

    Try to infer what the the underlying nucleotide change explains the observed AA change. This method exploits
    the heuristic that nearly all residue changes are due to SNPs (Pines and Winkler, 2016). However, this is not
    always the case, so these inferences are named inferred_nuc_snps instead of nuc_snps in the database.

    :param nt_sequence: sequence for the gene of interest
    :param gene_start: start position (should be zero indexed)
    :param wt_residue: residue at codon position
    :param position: residue position
    :param alt_aa: new residue
    :param aa_to_codons: mapping AA to codons
    :param codon_to_aa: codons to AA
    :return:
    """

    # get nucleotide window under consideration
    nucleotide_base_start = position * 3
    nucleotide_base_end = nucleotide_base_start + 3

    # this is the WT codon
    if position < 0 or position > len(nt_sequence):
        # just in case
        # the > condition is more likely
        raise AssertionError('Residue => nucleotide mapping is invalid!')

    wt_codon = nt_sequence[nucleotide_base_start:nucleotide_base_end]

    if position != 0:
        # if not the start codon, check to make sure the residues match
        # the start codon is a bit special in that E. coli will always use M instead of whatever codon is there
        # although a non-ATG codon will have much lower initiation efficiency
        if wt_residue != codon_to_aa[wt_codon]:
            raise AssertionError('Expected inferred and nominal WT residues to match, got: %s vs. %s'
                                 % (wt_residue, codon_to_aa[wt_codon]))

    # okay, now that we are here, we just need to get the codons for the alt amino acids, then see what mutations can
    # explain the residue switch

    if alt_aa not in aa_to_codons:
        # check for DB errors
        raise AssertionError('Unknown AA? %s' % alt_aa)

    hamming_dists = []

    for alt_codon in aa_to_codons[alt_aa]:
        hamming_dists.append((wt_codon, alt_codon, hamming_distance(wt_codon, alt_codon)))

    # we want to filter out the exact same codon if the mutation is synonymous
    hamming_dists = sorted(filter(lambda y: y[-1] > 0, hamming_dists), key=lambda x: x[-1])
    min_dist_counter = defaultdict(int)

    if len(hamming_dists) > 0:
        # if there are any surviving comparison tuples
        minimum_distance = hamming_dists[0][-1]
    else:
        # I don't think this should ever happen, but this will trigger a NullPointerException below if I am wrong
        minimum_distance = None

    for x, y, d in hamming_dists:
        if d == minimum_distance:
            min_dist_counter[d] += 1

    # only a single occurrence
    if min_dist_counter.get(minimum_distance, 0) == 1:
        # now we need to generate the 'inferred_nuc_snps' mutations
        position = nucleotide_base_start
        (wt_codon, alt_codon, _) = hamming_dists[0]
        output = []
        for wt_codon_base, alt_codon_base in zip(wt_codon, alt_codon):
            if wt_codon_base != alt_codon_base:
                # these codon positions do not match
                # generate mutations from it
                # assume gene_start is zero-indexed already
                output.append((position + gene_start, wt_codon_base, alt_codon_base))
            position += 1

        # and we're done.
        return output
    else:
        return []


def infer_residue_nucleotide_changes(cursor):
    """

    Handles the extraction of all AA residue changes and the inference of the underlying nucleotide changes.
    Returns a list of tuples for insertion into the resistome.annotations table.

    :param cursor:
    :return:
    """

    codon_to_aa, aa_to_codon = sql_interface.get_genetic_code()

    cursor.execute('SELECT resistome.annotations.annotation, '
                   'resistome.annotations.gene_id, '
                   'gene_locations.start,'
                   'nt_sequences.nt_seq '
                   'FROM resistome.annotations '
                   'INNER JOIN resistome.mutations on resistome.annotations.gene_id = resistome.mutations.gene_id '
                   'INNER JOIN genes ON genes.accession = resistome.mutations.name '
                   'INNER JOIN gene_locations ON gene_locations.gene_id = genes.gene_id '
                   'INNER JOIN nt_sequences on nt_sequences.gene_id = genes.gene_id '
                   'WHERE resistome.annotations.mutation_type = %s '
                   'AND resistome.annotations.valid IS TRUE', ('aa_snps',))

    grouped_updates = defaultdict(list)
    counter = 0
    for record in cursor:
        for aa_mutation in record['annotation']['aa_snps']:
            wt_aa = aa_mutation[1]
            position = aa_mutation[0]
            alt_aa = aa_mutation[-1]

            gene_id = record['gene_id']

            inferred_nuc_snps = infer_minimal_snps_from_aa(nt_sequence=record['nt_seq'],
                                                           gene_start=record['start'],
                                                           wt_residue=wt_aa,
                                                           position=position,
                                                           alt_aa=alt_aa,
                                                           aa_to_codons=aa_to_codon,
                                                           codon_to_aa=codon_to_aa)

            grouped_updates[gene_id].extend(inferred_nuc_snps)

    import json

    tuples_to_insert = []
    for gene_id, inferred_tuple_list in grouped_updates.items():
        # valid by definition.
        if len(inferred_tuple_list) == 0:
            continue

        tuples_to_insert.append((gene_id, 'inferred_nuc_snps', json.dumps({'inferred_nuc_snps': inferred_tuple_list}),
                                 True, None))

    return tuples_to_insert


if __name__ == '__main__':
    pass
