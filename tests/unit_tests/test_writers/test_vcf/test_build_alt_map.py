import os

import pandas as pd
import pytest

from CRISPResso2.CRISPRessoCOREResources import find_indels_substitutions
from CRISPResso2.writers import vcf

# ----------------------------- helpers -----------------------------

# This seq is altered for every test case; tests should not depend on its content.
REF_SEQ = "ATGCGTACGATCGTACGTAGCTAGCTAGCGTAGCTAGCTA"  # 40 bp
REF_POSITIONS = list(range(len(REF_SEQ)))


def create_df_alleles(*refs_alns):
    payloads = []
    for ref, *aln in refs_alns:
        if isinstance(aln, list):
            if len(aln) == 1:
                aln = aln[0]
                num_reads = 1
                ref_name = 'Reference'
            elif len(aln) == 2:
                aln, num_reads = aln
                ref_name = 'Reference'
            elif len(aln) == 3:
                aln, num_reads, ref_name = aln
        payload = find_indels_substitutions(aln, ref, list(range(len(ref)))).__dict__
        payload['Reference_Sequence'] = ref
        payload['Aligned_Sequence'] = aln
        payload['#Reads'] = num_reads
        payload['Reference_Name'] = ref_name
        payloads += [payload]

    for payload in payloads:
        payload['n_inserted'] = payload['insertion_n']
        payload['n_deleted'] = payload['deletion_n']
        payload['n_mutated'] = payload['substitution_n']

    return pd.DataFrame(payloads)


# ----------------------------- alteration functions -----------------------------
# These functions are used throughout the tests to alter the base REF_SEQ before calling build_edit_counts.

def make_unmodified(reads=1, ref_name="Reference"):
    """Unmodified read identical to REF_SEQ."""
    return (REF_SEQ, REF_SEQ, reads, ref_name)


def make_sub(position, alt_base, reads=1, ref_name="Reference"):
    """Single‑nucleotide substitution at 0‑based ref index 'position'."""
    aligned = REF_SEQ[:position] + alt_base + REF_SEQ[position + 1 :]
    return (REF_SEQ, aligned, reads, ref_name)


def make_del(start, end, reads=1, ref_name="Reference"):
    """Deletion [start, end) in reference coordinates (end exclusive)."""
    deletion_len = end - start
    aligned = REF_SEQ[:start] + "-"*deletion_len + REF_SEQ[end:]  # aligned sequence without the deleted block
    return (REF_SEQ, aligned, reads, ref_name)


def make_ins(after_index, ins_seq, reads=1, ref_name="Reference"):
    """
    Insertion of 'ins_seq' BETWEEN after_index and after_index+1 (0‑based).
    For build_edit_counts, insertion_coordinates use the right‑anchor ref index (after_index+1),
    and the aligned start is exactly after_index+1 when there are no prior edits in the row.
    """
    aligned_start = after_index + 1
    aligned = REF_SEQ[:aligned_start] + ins_seq + REF_SEQ[aligned_start:]
    updated_ref_seq = REF_SEQ[:aligned_start] + "-"*len(ins_seq) + REF_SEQ[aligned_start:]
    return (updated_ref_seq, aligned, reads, ref_name)


# ----------------------------- core mixed case -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # Mixed: SNP at index 9, 1‑bp del at index 19, 10‑bp del at index 19, 2‑bp insertion after 30, plus an unmodified that must be skipped.
        (
            [
                make_sub(9, "G", reads=1),
                make_del(19, 20, reads=1),
                make_ins(30, "GG", reads=1),
                make_del(19, 29, reads=1),
                make_unmodified(reads=2),
            ],
            {"Reference": (1, 1)},
            {
                (1, 10, REF_SEQ[9], "G"): 1,
                (1, 19, REF_SEQ[18:20], REF_SEQ[18]): 1,       # 1‑bp del
                (1, 19, REF_SEQ[18:29], REF_SEQ[18]): 1,       # 10‑bp del
                (1, 31, REF_SEQ[30], REF_SEQ[30] + "GG"): 1,   # insertion
            },
        ),
    ],
    ids=["mixed_edits_happy_path"],
)
def test_build_edit_counts_mixed(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


# ----------------------------- substitutions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single SNP
        (
            [make_sub(9, "G", reads=1)],
            {"Reference": (1, 1)},
            {(1, 10, REF_SEQ[9], "G"): 1},
        ),
        # merge same SNP (reads sum)
        (
            [make_sub(9, "G", reads=2), make_sub(9, "G", reads=3)],
            {"Reference": (1, 1)},
            {(1, 10, REF_SEQ[9], "G"): 5},
        ),
        # split different alt bases
        (
            [make_sub(9, "G", reads=1), make_sub(9, "T", reads=2)],
            {"Reference": (1, 1)},
            {(1, 10, REF_SEQ[9], "G"): 1, (1, 10, REF_SEQ[9], "T"): 2},
        ),
    ],
    ids=["sub_single", "sub_merge_same_base", "sub_split_diff_base"],
)
def test_build_edit_counts_substitutions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


# ----------------------------- deletions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single 1‑bp deletion at start=19
        (
            [make_del(19, 20, reads=1)],
            {"Reference": (1, 1)},
            {(1, 19, REF_SEQ[18:20], REF_SEQ[18]): 1},
        ),
        # merge same‑length deletion at same key (reads sum)
        (
            [make_del(19, 20, reads=2), make_del(19, 20, reads=3)],
            {"Reference": (1, 1)},
            {(1, 19, REF_SEQ[18:20], REF_SEQ[18]): 5},
        ),
        # two different lengths at same position -> two entries
        (
            [make_del(19, 20, reads=1), make_del(19, 29, reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 19, REF_SEQ[18:20], REF_SEQ[18]): 1,
                (1, 19, REF_SEQ[18:29], REF_SEQ[18]): 1,
            },
        ),
        # deletion at second element
        (
            [make_del(1, 2, reads=1)],
            {'Reference': (1, 1)},
            {(1, 1, REF_SEQ[0:2], REF_SEQ[0]): 1},
        ),
    ],
    ids=["del_single", "del_merge_same_len", "del_two_lengths_same_key", "del_second_element"],
)
def test_build_edit_counts_deletions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        (
            [make_del(39, 40, reads=1)],
            {'Reference': (1, 1)},
            {(1, 39, REF_SEQ[38:40], REF_SEQ[38]): 1},
        ),
        (
            [make_del(38, 40, reads=1)],
            {"Reference": (1, 1)},
            {(1, 38, REF_SEQ[37:40], REF_SEQ[37]): 1},
        ),
    ],
    ids=['last_element', 'second_to_last_element']
)
def test_build_edit_counts_deletion_at_end(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


def test_build_edit_counts_deletion_start_at_zero_should_anchor_correctly():
    """When a deletion occurs at the start of a sequence, REF includes the base after the deletion.

    Source: https://bioinformatics.stackexchange.com/questions/2476/how-to-represent-a-deletion-at-position-1-in-a-vcf-file
    """
    df = create_df_alleles(make_del(0, 3, reads=1))  # delete first 3 bases
    amplicon_positions = {"Reference": (1, 1)}
    out = vcf.build_edit_counts(df, amplicon_positions)
    # REF = deleted + following_base, ALT = following_base
    assert out == {(1, 1, REF_SEQ[0:4], REF_SEQ[3]): 1}


def test_build_edit_counts_deletion_start_at_zero_pos_offset_single_base():
    """Single-base deletion at start=0 with pos=100 should give VCF POS=100, not 99."""
    ref = 'GATTACA'
    aln = '-ATTACA'
    df = create_df_alleles((ref, aln))
    amplicon_positions = {"Reference": ("chr1", 100)}
    out = vcf.build_edit_counts(df, amplicon_positions)
    # REF = deleted + following_base, ALT = following_base
    assert out == {("chr1", 100, "GA", "A"): 1}


def test_build_edit_counts_deletion_start_at_zero_pos_offset_multi_base():
    """Multi-base deletion at start=0 with pos=100 should give VCF POS=100, not 99."""
    df = create_df_alleles(make_del(0, 3, reads=1))  # delete first 3 bases
    amplicon_positions = {"Reference": ("chr1", 100)}
    out = vcf.build_edit_counts(df, amplicon_positions)
    # REF = deleted + following_base, ALT = following_base
    assert out == {("chr1", 100, REF_SEQ[0:4], REF_SEQ[3]): 1}


def test_build_edit_counts_deletion_start():
    ref1 = 'GATTACA'
    aln1 = '-ATTACA'

    df = create_df_alleles((ref1, aln1))
    amplicon_positions = {"Reference": (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    assert edit_counts == {(1, 1, 'GA', 'A'): 1}


def test_build_edit_counts_multi_deletion_start():
    ref = 'AACCTTGG'
    df = create_df_alleles(
        (ref, '-ACCTTGG'),
        (ref, '--CCTTGG'),
        (ref, '---CTTGG'),
        (ref, '----TTGG'),
        (ref, '-----TGG'),
        (ref, '------GG'),
        (ref, '-------G'),
    )

    amplicon_positions = {"Reference": (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    assert edit_counts == {
        (1, 1, 'AA', 'A'): 1,
        (1, 1, 'AAC', 'C'): 1,
        (1, 1, 'AACC', 'C'): 1,
        (1, 1, 'AACCT', 'T'): 1,
        (1, 1, 'AACCTT', 'T'): 1,
        (1, 1, 'AACCTTG', 'G'): 1,
        (1, 1, 'AACCTTGG', 'G'): 1,
    }

    # Validate VCF output: 7 biallelic lines
    temp_vcf_path = 'multi_deletion_start.vcf'
    num_vcf_rows = vcf.write_vcf_from_edits(edit_counts, len(df), {}, temp_vcf_path)
    assert num_vcf_rows == 7

    with open(temp_vcf_path) as fh:
        vcf_contents = fh.read()
    assert 'AA\tA' in vcf_contents
    assert 'AAC\tC' in vcf_contents
    assert 'AACC\tC' in vcf_contents
    assert 'AACCT\tT' in vcf_contents
    assert 'AACCTT\tT' in vcf_contents
    assert 'AACCTTG\tG' in vcf_contents
    assert 'AACCTTGG\tG' in vcf_contents

    os.remove(temp_vcf_path)


def test_build_edit_counts_multi_deletion_start_and_middle():
    ref = 'AACCTTGG'
    df = create_df_alleles(
        (ref, '-ACCTTGG'),
        (ref, '--CCTTGG'),
        (ref, 'A--CTTGG'),   # middle deletion (2 reads)
        (ref, 'A--CTTGG'),
        (ref, '----TTGG'),
        (ref, '-----TGG'),
        (ref, '------GG'),
        (ref, '-------G'),   # (2 reads)
        (ref, '-------G'),
    )

    amplicon_positions = {"Reference": (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    assert edit_counts == {
        (1, 1, 'AA', 'A'): 1,          # del_start 1bp
        (1, 1, 'AAC', 'C'): 1,         # del_start 2bp
        (1, 1, 'AAC', 'A'): 2,         # middle del 'AC' (2 reads)
        (1, 1, 'AACCT', 'T'): 1,       # del_start 4bp
        (1, 1, 'AACCTT', 'T'): 1,      # del_start 5bp
        (1, 1, 'AACCTTG', 'G'): 1,     # del_start 6bp
        (1, 1, 'AACCTTGG', 'G'): 2,    # del_start 7bp (2 reads)
    }

    # Validate VCF output
    temp_vcf_path = 'multi_deletion_start_and_middle.vcf'
    num_vcf_rows = vcf.write_vcf_from_edits(edit_counts, len(df), {}, temp_vcf_path)
    assert num_vcf_rows == 7

    with open(temp_vcf_path) as fh:
        vcf_contents = fh.read()
    assert 'AA\tA' in vcf_contents
    assert 'AAC\tC' in vcf_contents
    assert 'AAC\tA' in vcf_contents
    assert 'AACCT\tT' in vcf_contents
    assert 'AACCTT\tT' in vcf_contents
    assert 'AACCTTG\tG' in vcf_contents
    assert 'AACCTTGG\tG' in vcf_contents

    os.remove(temp_vcf_path)


# ----------------------------- insertions -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        # single insertion of "GG" after index 30
        (
            [make_ins(30, "GG", reads=1)],
            {"Reference": (1, 1)},
            {(1, 31, REF_SEQ[30], REF_SEQ[30] + "GG"): 1},
        ),
        # merge same inserted string and key (reads sum)
        (
            [make_ins(30, "GG", reads=2), make_ins(30, "GG", reads=3)],
            {"Reference": (1, 1)},
            {(1, 31, REF_SEQ[30], REF_SEQ[30] + "GG"): 5},
        ),
        # split different inserted strings at same anchor (before normalization)
        # "GG" at anchor 30(T): ins[-1]='G' != 'T' -> stays at pos 31
        # "T" at anchor 30(T): ins[-1]='T' == 'T' -> left-normalizes to anchor 29(G), pos 30
        (
            [make_ins(30, "GG", reads=1), make_ins(30, "T", reads=1)],
            {"Reference": (1, 1)},
            {
                (1, 31, REF_SEQ[30], REF_SEQ[30] + "GG"): 1,
                (1, 30, REF_SEQ[29], REF_SEQ[29] + "T"): 1,
            },
        ),
    ],
    ids=["ins_single", "ins_merge_same_seq", "ins_split_diff_seq"],
)
def test_build_edit_counts_insertions(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


def test_build_edit_counts_two_insertions_same_read():
    """Two insertions in the same read should both extract the correct bases.

    Amplicon: ATGCATGC (8bp)
    Read has CC inserted after pos 1, and GG inserted after pos 3.
    Aligned ref: AT--GC--ATGC
    Aligned read: ATCCGCGGATGC
    """
    ref = 'AT--GC--ATGC'
    aln = 'ATCCGCGGATGC'
    df = create_df_alleles((ref, aln, 1))
    amplicon_positions = {"Reference": ("chr1", 1)}
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == {
        ("chr1", 2, "T", "TCC"): 1,   # insertion of CC, anchor at ref pos 1 (T)
        ("chr1", 4, "C", "CGG"): 1,   # insertion of GG, anchor at ref pos 3 (C)
    }


# ----------------------------- multiple amplicons & coordinate offsets -----------------------------
@pytest.mark.parametrize(
    "rows, amplicon_positions, expected",
    [
        (
            # One SNP on amplicon A, one insertion on amplicon B, ensure independence and absolute coordinates.
            [make_sub(2, "T", reads=1, ref_name="A"),
             make_ins(0, "G", reads=2, ref_name="B")],
            {"A": (1, 1), "B": (2, 1000)},
            {
                (1, 3, REF_SEQ[2], "T"): 1,
                (2, 1000, REF_SEQ[0], REF_SEQ[0] + "G"): 2,
            },
        ),
    ],
    ids=["two_amplicons_offsets"],
)
def test_build_edit_counts_multi_amplicon_and_offsets(rows, amplicon_positions, expected):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


# ----------------------------- unmodified rows are skipped -----------------------------

@pytest.mark.parametrize(
    "rows, amplicon_positions, expected_size",
    [
        ([make_unmodified(reads=5)], {"Reference": (1, 1)}, 0),
        ([make_unmodified(), make_sub(0, "C", 1)], {"Reference": (1, 1)}, 1),
    ],
    ids=["only_unmodified", "mixed_with_unmodified"],
)
def test_build_edit_counts_skips_unmodified(rows, amplicon_positions, expected_size):
    df = create_df_alleles(*rows)
    out = vcf.build_edit_counts(df, amplicon_positions)
    assert len(out) == expected_size

# ----------------------------- fidelity test with real-like alignment data -----------------------------

def test_build_edit_counts_fidelity_like_real_row():
    ref_seq = "ATGCGTACGA---TCGTACGTAGCTAGCTAGCGTAGCTAGCTA"
    aln_seq = "ATGCGTACGAGGGTCGTACGTAGCTAGCTAGCGTAGCTAGCTA"
    num_reads = 7
    ins_right_anchor = 9

    df = create_df_alleles((ref_seq, aln_seq, num_reads))
    amplicon_positions = {"Reference": (1, 500)}

    anchor = ref_seq[ins_right_anchor]
    expected = {
        (1, 500 + ins_right_anchor, anchor, anchor + "GGG"): num_reads,
    }

    out = vcf.build_edit_counts(df, amplicon_positions)
    assert out == expected


def test_aln_to_edit_counts_ins_del_same_pos():
    ref1 = 'AATGCGTAC'
    aln1 = 'AA--CGTAC'

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AATTTGCGTAC'

    df = create_df_alleles((ref1, aln1), (ref2, aln2))
    amplicon_positions = {"Reference": (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    assert edit_counts == {
        (1, 2, 'ATG', 'A'): 1,    # deletion of TG
        (1, 2, 'A', 'ATT'): 1,    # insertion of TT
    }

    temp_vcf_path = 'aln_to_edit_counts_ins_del_same_pos.vcf'
    num_reads = 5
    num_vcf_rows = vcf.write_vcf_from_edits(edit_counts, num_reads, {'Reference': 9}, temp_vcf_path)
    assert num_vcf_rows == 2

    with open(temp_vcf_path) as fh:
        vcf_contents = fh.read()
    assert 'ATG\tA' in vcf_contents
    assert 'A\tATT' in vcf_contents

    os.remove(temp_vcf_path)


def test_aln_to_edit_counts_ins_then_del():
    ref1 = 'AATTT---GCGTAC'
    aln1 = 'AATTTCCCGCG---'

    df = create_df_alleles((ref1, aln1))
    amplicon_positions = {'Reference': (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    assert edit_counts == {
        (1, 5, 'T', 'TCCC'): 1,     # insertion of CCC
        (1, 8, 'GTAC', 'G'): 1,     # deletion of TAC
    }


def test_aln_to_edit_counts_to_vcf():
    ref1 = 'AATGCGTAC'
    aln1 = 'AATGCG-AC'
    #             ^ Interested in this deletion across each of these examples

    ref2 = 'AA--TGCGTAC'
    aln2 = 'AAGGTGCG-AC'

    ref3 = 'AATGCGTAC'
    aln3 = 'AA-GCG-AC'

    ref4 = 'AATG--CGTAC'
    aln4 = 'AA-GAACG-AC'

    ref5 = 'AA--TGCGTAC'
    aln5 = 'AAGGTGCG-GG'

    ref6 = 'AATGCGTAC'
    aln6 = '-------AC'

    df = create_df_alleles(
        (ref1, aln1),
        (ref2, aln2),
        (ref3, aln3),
        (ref4, aln4),
        (ref5, aln5),
        (ref6, aln6),
    )

    amplicon_positions = {"Reference": (1, 1)}
    edit_counts = vcf.build_edit_counts(df, amplicon_positions)

    # deletion at position 7 (1-based) in each example above, should occur 5 times
    assert edit_counts[(1, 6, 'GT', 'G')] == 5
    # insertion of GG occurs 2 times in 2, 5 and deletion of T occurs 2 times in 3, 4
    assert edit_counts[(1, 2, 'A', 'AGG')] == 2
    assert edit_counts[(1, 2, 'AT', 'A')] == 2
    # insertion of AA occurs 1 time in 4
    assert edit_counts[(1, 4, 'G', 'GAA')] == 1
    # substitution of A -> G occurs 1 time in 5
    assert edit_counts[(1, 8, 'A', 'G')] == 1
    # substitution of C -> G occurs 1 time in 5
    assert edit_counts[(1, 9, 'C', 'G')] == 1
    # deletion at position 1 (1-based) occurs 1 time in 6
    assert edit_counts[(1, 1, 'AATGCGTA', 'A')] == 1

    temp_vcf_path = 'aln_to_edit_counts_to_vcf.vcf'
    num_reads = 5
    num_vcf_rows = vcf.write_vcf_from_edits(edit_counts, num_reads, {'Reference': 9}, temp_vcf_path)
    # Biallelic: 7 records (position 2 has 2 records, rest have 1 each)
    assert num_vcf_rows == 7

    with open(temp_vcf_path, 'r') as fh:
        vcf_contents = fh.read()

    assert '\t'.join(('1', '6', '.', 'GT', 'G', '.', 'PASS', f'AF={5 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '2', '.', 'AT', 'A', '.', 'PASS', f'AF={2 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '2', '.', 'A', 'AGG', '.', 'PASS', f'AF={2 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '4', '.', 'G', 'GAA', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '8', '.', 'A', 'G', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '9', '.', 'C', 'G', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents
    assert '\t'.join(('1', '1', '.', 'AATGCGTA', 'A', '.', 'PASS', f'AF={1 / num_reads:.3f}')) in vcf_contents

    os.remove(temp_vcf_path)
