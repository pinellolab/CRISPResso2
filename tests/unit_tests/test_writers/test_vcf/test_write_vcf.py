import os
from types import SimpleNamespace

import pytest

from CRISPResso2 import CRISPRessoShared
from CRISPResso2.writers import vcf


# ----------------------------- write_vcf_from_edits -----------------------------

def test_write_vcf_from_edits_header_and_ordering_no_samples():
    """Tests that write_vcf_from_edits produces the expected header and ordering of records."""
    edit_counts = {
        ("1", 20, "AG", "A"): 2,   # deletion
        ("1", 10, "A", "T"): 3,    # substitution
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_no_samples.vcf")
    count = vcf.write_vcf_from_edits(edit_counts, num_reads=10, amplicon_lens={}, vcf_path=temp_vcf_path)

    with open(temp_vcf_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines[0].startswith("##fileformat=VCFv4.")
    assert lines[1] == "##source=CRISPResso2"
    assert lines[2].startswith("##INFO=<ID=AF")
    assert lines[3].startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # Records in ascending chrom/pos order
    assert lines[4] == "1\t10\t.\tA\tT\t.\tPASS\tAF=0.300"
    assert lines[5] == "1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200"

    assert len(lines) == 6
    assert count == 2
    os.remove(temp_vcf_path)


def test_write_vcf_from_edits_header_and_ordering_with_contigs():
    """Test that write_vcf_from_edits produces contig headers and the expected record."""
    edit_counts = {
        ("X", 2, "A", "AG"): 2,   # insertion
    }
    temp_vcf_path = os.path.join(os.path.dirname(__file__), "temp_test_with_contigs.vcf")
    count = vcf.write_vcf_from_edits(
        edit_counts, num_reads=10,
        amplicon_lens={"Ref1": 100, "Ref2": 200},
        vcf_path=temp_vcf_path,
    )

    with open(temp_vcf_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines[3] == "##contig=<ID=Ref1,length=100>"
    assert lines[4] == "##contig=<ID=Ref2,length=200>"
    assert lines[5] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    assert lines[6] == "X\t2\t.\tA\tAG\t.\tPASS\tAF=0.200"
    assert len(lines) == 7
    assert count == 1
    os.remove(temp_vcf_path)


def test_write_vcf_from_edits_bad_denominator_raises():
    for bad_num_reads in [0, None, -5]:
        with pytest.raises(CRISPRessoShared.BadParameterException):
            vcf.write_vcf_from_edits(
                {("1", 10, "A", "T"): 1},
                bad_num_reads,
                {},
                "/dev/null",
            )


# ----------------------------- write_vcf_file (orchestrator) -----------------------------

def _stub_edit_counts():
    return {
        ("1", 10, "A", "T"): 3,
        ("1", 20, "AG", "A"): 2,
    }

@pytest.mark.parametrize(
    "ref_names,ref_lens,coords_str,expected_record_count,expected_contig_lines",
    [
        (["RefA"], {"RefA": 100}, "1:1", 2, 1),
        (["RefA", "RefB"], {"RefA": 100, "RefB": 200}, "1:1,2:1", 2, 2),
    ],
)
def test_write_vcf_file_smoke(tmp_path, monkeypatch, ref_names, ref_lens, coords_str, expected_record_count, expected_contig_lines):
    monkeypatch.setattr(vcf, "build_edit_counts", lambda df, amplicon_positions: _stub_edit_counts())

    import pandas as pd
    df = pd.DataFrame([{"#Reads": 7}, {"#Reads": 3}])  # num_reads = 10

    args = SimpleNamespace(amplicon_coordinates=coords_str)
    out_path = tmp_path / "test.vcf"

    num_vcf_rows = vcf.write_vcf_file(df, ref_names, ref_lens, args, str(out_path))

    with open(out_path, "r") as f:
        lines = [line.strip() for line in f]

    assert num_vcf_rows == 2

    expected_header_lines = 4 + expected_contig_lines
    assert len(lines) == expected_header_lines + expected_record_count

    assert lines[expected_header_lines].startswith("1\t10\t.\tA\tT\t.\tPASS\tAF=0.300")
    assert lines[expected_header_lines + 1].startswith("1\t20\t.\tAG\tA\t.\tPASS\tAF=0.200")


# ----------------------------- _left_normalize_deletion -----------------------------


def _simple_ref(seq):
    """Return (ref_positions, ref_str) for an ungapped reference."""
    return list(range(len(seq))), seq


class TestLeftNormalizeDeletion:
    """Direct tests for vcf._left_normalize_deletion."""

    def test_no_shift_distinct_bases(self):
        """ATGC: deleting G(2,3) - base before is T, last deleted is G -> no shift."""
        rp, rs = _simple_ref("ATGC")
        assert vcf._left_normalize_deletion(2, 3, rp, rs) == (2, 3)

    def test_single_base_shift(self):
        """ATCC: deleting C at pos 3 -> shifts left to pos 2 (both C)."""
        #        0123
        # ref:   ATCC
        # base before pos 3 = C(pos 2), last deleted = C(pos 2 after shift) -> shift once
        # base before pos 2 = T(pos 1), last deleted = C(pos 2) -> stop
        rp, rs = _simple_ref("ATCC")
        assert vcf._left_normalize_deletion(3, 4, rp, rs) == (2, 3)

    def test_multi_base_shift(self):
        """AAAAATG: deleting A at pos 4 -> shifts all the way to pos 0."""
        #        0123456
        # ref:   AAAAATG
        rp, rs = _simple_ref("AAAAATG")
        assert vcf._left_normalize_deletion(4, 5, rp, rs) == (0, 1)

    def test_multi_base_deletion_shift(self):
        """ACACACG: deleting AC at (4,6) -> shifts to (0,2)."""
        #        0123456
        # ref:   ACACACG
        # del (4,6) = AC. base before=C(3), last del=C(5) -> equal, shift to (3,5)
        # del (3,5) = AC. base before=A(2), last del=A(4) -> equal, shift to (2,4)
        # del (2,4) = AC. base before=C(1), last del=C(3) -> equal, shift to (1,3)
        # del (1,3) = CA. base before=A(0), last del=A(2) -> equal, shift to (0,2)
        # start=0, stop
        rp, rs = _simple_ref("ACACACG")
        assert vcf._left_normalize_deletion(4, 6, rp, rs) == (0, 2)

    def test_already_at_position_zero(self):
        """Deletion already at start of sequence - can't shift further."""
        rp, rs = _simple_ref("ATGC")
        assert vcf._left_normalize_deletion(0, 1, rp, rs) == (0, 1)

    def test_homopolymer_shifts_to_start(self):
        """AAAAT: deleting A at pos 3 -> shifts to pos 0."""
        rp, rs = _simple_ref("AAAAT")
        assert vcf._left_normalize_deletion(3, 4, rp, rs) == (0, 1)

    def test_partial_homopolymer_shift(self):
        """TAAAC: deleting A at pos 3 -> shifts to pos 1 (T blocks further)."""
        rp, rs = _simple_ref("TAAAC")
        assert vcf._left_normalize_deletion(3, 4, rp, rs) == (1, 2)

    def test_dinucleotide_repeat_multi_del(self):
        """GATGATGATC: deleting GAT at (6,9) -> shifts to (0,3)."""
        #        0123456789
        # ref:   GATGATGATC
        rp, rs = _simple_ref("GATGATGATC")
        assert vcf._left_normalize_deletion(6, 9, rp, rs) == (0, 3)

    def test_no_shift_when_boundary_bases_differ(self):
        """ACGT: deleting CG at (1,3) - base before=A(0), last del=G(2) -> no shift."""
        rp, rs = _simple_ref("ACGT")
        assert vcf._left_normalize_deletion(1, 3, rp, rs) == (1, 3)


# ----------------------------- _left_normalize_insertion -----------------------------


class TestLeftNormalizeInsertion:
    """Direct tests for vcf._left_normalize_insertion."""

    def test_no_shift_distinct_bases(self):
        """ATGC: inserting 'A' at anchor pos 2 (G). Last ins base A != G -> no shift."""
        rp, rs = _simple_ref("ATGC")
        anchor, ins = vcf._left_normalize_insertion(2, "A", rp, rs)
        assert anchor == 2
        assert ins == "A"

    def test_single_base_shift_homopolymer(self):
        """ATCC: inserting 'C' at anchor pos 2 (C). Shifts left to pos 1 (T != C -> stop)."""
        #        0123
        # ref:   ATCC
        # anchor=2(C), ins='C'. ins[-1]='C' == ref[2]='C' -> ins='C'+'', anchor=1
        # anchor=1(T), ins='C'. ins[-1]='C' != ref[1]='T' -> stop
        rp, rs = _simple_ref("ATCC")
        anchor, ins = vcf._left_normalize_insertion(2, "C", rp, rs)
        assert anchor == 1
        assert ins == "C"

    def test_homopolymer_shifts_to_zero(self):
        """ACCCG: inserting 'C' at anchor pos 3. Shifts to anchor 0 (A != C -> stop)."""
        #        01234
        # ref:   ACCCG
        rp, rs = _simple_ref("ACCCG")
        anchor, ins = vcf._left_normalize_insertion(3, "C", rp, rs)
        assert anchor == 0
        assert ins == "C"

    def test_already_at_position_zero(self):
        """ACGT: inserting 'A' at anchor pos 0. Can't shift further."""
        rp, rs = _simple_ref("ACGT")
        anchor, ins = vcf._left_normalize_insertion(0, "A", rp, rs)
        assert anchor == 0
        assert ins == "A"

    def test_multi_base_insertion_dinucleotide_repeat(self):
        """ACACG: inserting 'AC' at anchor pos 3 (C). Shifts to anchor 0.

        anchor=3(C), ins='AC'. ins[-1]='C' == ref[3]='C' -> ins='C'+'A'='CA', anchor=2
        anchor=2(A), ins='CA'. ins[-1]='A' == ref[2]='A' -> ins='A'+'C'='AC', anchor=1
        anchor=1(C), ins='AC'. ins[-1]='C' == ref[1]='C' -> ins='C'+'A'='CA', anchor=0
        anchor=0, stop.
        """
        rp, rs = _simple_ref("ACACG")
        anchor, ins = vcf._left_normalize_insertion(3, "AC", rp, rs)
        assert anchor == 0
        assert ins == "CA"

    def test_no_shift_when_last_ins_base_differs(self):
        """ACGT: inserting 'T' at anchor pos 2 (G). T != G -> no shift."""
        rp, rs = _simple_ref("ACGT")
        anchor, ins = vcf._left_normalize_insertion(2, "T", rp, rs)
        assert anchor == 2
        assert ins == "T"

    def test_partial_shift_multi_base(self):
        """TAAAC: inserting 'A' at anchor pos 3 (A). Shifts to anchor 1 (T blocks).

        anchor=3(A), ins='A'. ins[-1]='A' == ref[3]='A' -> ins='A', anchor=2
        anchor=2(A), ins='A'. ins[-1]='A' == ref[2]='A' -> ins='A', anchor=1
        anchor=1(A), ins='A'. ins[-1]='A' == ref[1]='A' -> ins='A', anchor=0
        anchor=0(T), ins='A'. ins[-1]='A' != ref[0]='T' -> stop.
        """
        rp, rs = _simple_ref("TAAAC")
        anchor, ins = vcf._left_normalize_insertion(3, "A", rp, rs)
        assert anchor == 0
        assert ins == "A"

    def test_trinucleotide_repeat(self):
        """GATGATGATC: inserting 'GAT' at anchor pos 8 (T). Shifts to anchor 0.

        anchor=8(T), ins='GAT'. ins[-1]='T' == ref[8]='T' -> ins='T'+'GA'='TGA', anchor=7
        anchor=7(A), ins='TGA'. ins[-1]='A' == ref[7]='A' -> ins='A'+'TG'='ATG', anchor=6
        anchor=6(G), ins='ATG'. ins[-1]='G' == ref[6]='G' -> ins='G'+'AT'='GAT', anchor=5
        anchor=5(T), ins='GAT'. ins[-1]='T' == ref[5]='T' -> ins='T'+'GA'='TGA', anchor=4
        anchor=4(A), ins='TGA'. ins[-1]='A' == ref[4]='A' -> ins='A'+'TG'='ATG', anchor=3
        anchor=3(G), ins='ATG'. ins[-1]='G' == ref[3]='G' -> ins='G'+'AT'='GAT', anchor=2
        anchor=2(T), ins='GAT'. ins[-1]='T' == ref[2]='T' -> ins='T'+'GA'='TGA', anchor=1
        anchor=1(A), ins='TGA'. ins[-1]='A' == ref[1]='A' -> ins='A'+'TG'='ATG', anchor=0
        anchor=0, stop.
        """
        rp, rs = _simple_ref("GATGATGATC")
        anchor, ins = vcf._left_normalize_insertion(8, "GAT", rp, rs)
        assert anchor == 0
        assert ins == "ATG"


# ----------------------------- _edits_from_insertions left-normalization ---------------------


def _make_insertion_row(amplicon, anchor_ref_pos, ins_bases, reads=1):
    """Build a minimal allele-row dict for an insertion-only read.

    Constructs a synthetic aligned sequence with the insertion spliced in,
    and computes ref_positions with -1 entries for inserted bases.
    """
    ref_positions = []
    aligned_ref = []
    aligned_read = []

    for ref_pos in range(len(amplicon)):
        ref_positions.append(ref_pos)
        aligned_ref.append(amplicon[ref_pos])
        aligned_read.append(amplicon[ref_pos])

        if ref_pos == anchor_ref_pos:
            aligned_start = len(aligned_read)
            for base in ins_bases:
                ref_positions.append(-1)
                aligned_ref.append("-")
                aligned_read.append(base)

    return {
        "ref_positions": ref_positions,
        "Reference_Sequence": "".join(aligned_ref),
        "Aligned_Sequence": "".join(aligned_read),
        "#Reads": reads,
        "n_deleted": 0,
        "n_inserted": len(ins_bases),
        "n_mutated": 0,
        "deletion_coordinates": [],
        "deletion_sizes": [],
        "Reference_Name": "Reference",
        "insertion_coordinates": [(anchor_ref_pos, aligned_start)],
        "insertion_sizes": [len(ins_bases)],
        "substitution_positions": [],
    }


def test_insertion_left_normalization_homopolymer():
    """Inserting 'C' at anchor pos 3 in ACCCG should left-normalize to anchor pos 0.

    VCF before normalization: POS=4, REF=C, ALT=CC
    VCF after normalization:  POS=1, REF=A, ALT=AC
    """
    row = _make_insertion_row("ACCCG", 3, "C")
    edits = list(vcf._edits_from_insertions(row, "chr1", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert chrom == "chr1"
    assert vcf_pos == 1, f"Expected left-normalized VCF pos 1, got {vcf_pos}"
    assert ref == "A", f"Expected REF=A, got {ref}"
    assert alt == "AC", f"Expected ALT=AC, got {alt}"


def test_insertion_no_normalization_needed():
    """Inserting 'T' at anchor pos 2 in ACGT: no shift (T != G)."""
    row = _make_insertion_row("ACGT", 2, "T")
    edits = list(vcf._edits_from_insertions(row, "chr1", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert vcf_pos == 3
    assert ref == "G"
    assert alt == "GT"


def test_insertion_left_normalization_dinucleotide():
    """Inserting 'AC' at anchor pos 3 in ACACG should left-normalize to anchor pos 0.

    VCF before: POS=4, REF=C, ALT=CAC
    VCF after:  POS=1, REF=A, ALT=ACA
    """
    row = _make_insertion_row("ACACG", 3, "AC")
    edits = list(vcf._edits_from_insertions(row, "chr1", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert vcf_pos == 1, f"Expected left-normalized VCF pos 1, got {vcf_pos}"
    assert ref == "A", f"Expected REF=A, got {ref}"
    assert alt == "ACA", f"Expected ALT=ACA, got {alt}"


# ----------------------------- _edits_from_insertions multi-insertion --------------------------


def _make_two_insertion_row(amplicon, after1, ins1, after2, ins2, reads=1):
    """Build a minimal allele-row dict for a read with two insertions.

    Insertions are placed after 0-based ref positions after1 and after2.
    after1 must be < after2.
    """
    ref_positions = []
    aligned_ref = []
    aligned_read = []
    insertion_coordinates = []
    insertion_sizes = []

    for ref_pos in range(len(amplicon)):
        ref_positions.append(ref_pos)
        aligned_ref.append(amplicon[ref_pos])
        aligned_read.append(amplicon[ref_pos])

        if ref_pos == after1:
            for base in ins1:
                ref_positions.append(-1)
                aligned_ref.append("-")
                aligned_read.append(base)

        if ref_pos == after2:
            for base in ins2:
                ref_positions.append(-1)
                aligned_ref.append("-")
                aligned_read.append(base)

    # Derive insertion_coordinates using find_indels_substitutions
    from CRISPResso2.CRISPRessoCOREResources import find_indels_substitutions
    r = find_indels_substitutions(
        "".join(aligned_read), "".join(aligned_ref),
        list(range(len(aligned_ref))),
    )

    return {
        "ref_positions": ref_positions,
        "Reference_Sequence": "".join(aligned_ref),
        "Aligned_Sequence": "".join(aligned_read),
        "#Reads": reads,
        "n_deleted": 0,
        "n_inserted": len(ins1) + len(ins2),
        "n_mutated": 0,
        "deletion_coordinates": [],
        "deletion_sizes": [],
        "Reference_Name": "Reference",
        "insertion_coordinates": r.insertion_coordinates,
        "insertion_sizes": r.insertion_sizes,
        "substitution_positions": [],
    }


def test_two_insertions_same_read():
    """Two insertions in one read should both extract the correct bases.

    Amplicon: ATGCATGC
    Insert CC after pos 1, GG after pos 3.
    Aligned ref:  AT--GC--ATGC
    Aligned read: ATCCGCGGATGC

    First insertion (CC): anchor=T(pos1), VCF POS=2, REF=T, ALT=TCC
    Second insertion (GG): anchor=C(pos3), VCF POS=4, REF=C, ALT=CGG
    """
    row = _make_two_insertion_row("ATGCATGC", 1, "CC", 3, "GG")
    edits = list(vcf._edits_from_insertions(row, "chr1", 1))
    assert len(edits) == 2

    chrom1, pos1, ref1, alt1, reads1 = edits[0]
    assert pos1 == 2, f"Expected VCF pos 2, got {pos1}"
    assert ref1 == "T", f"Expected REF=T, got {ref1}"
    assert alt1 == "TCC", f"Expected ALT=TCC, got {alt1}"

    chrom2, pos2, ref2, alt2, reads2 = edits[1]
    assert pos2 == 4, f"Expected VCF pos 4, got {pos2}"
    assert ref2 == "C", f"Expected REF=C (anchor at ref pos 3), got {ref2}"
    assert alt2 == "CGG", f"Expected ALT=CGG, got {alt2}"


# ----------------------------- _edits_from_deletions left-normalization ----------------------

# The FANC amplicon used in integration tests
_FANC_AMPLICON = (
    "CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATG"
    "GAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGG"
    "ACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG"
)


def _make_deletion_row(amplicon, del_start, del_end, reads=1):
    """Build a minimal allele-row dict for a deletion-only read."""
    ref_positions = list(range(len(amplicon)))
    return {
        "ref_positions": ref_positions,
        "Reference_Sequence": amplicon,
        "Aligned_Sequence": amplicon,
        "#Reads": reads,
        "n_deleted": del_end - del_start,
        "n_inserted": 0,
        "n_mutated": 0,
        "deletion_coordinates": [(del_start, del_end)],
        "deletion_sizes": [del_end - del_start],
        "Reference_Name": "Reference",
        "insertion_coordinates": [],
        "insertion_sizes": [],
        "substitution_positions": [],
    }


def test_deletion_left_normalization_two_base_shift():
    """Deletion (92,94) in FANC amplicon should left-normalize to VCF pos 91, REF=GCA, ALT=G.

    The aligner may right-shift deletions in repetitive regions (here ...CACCT...).
    Deleting bases 92-93 (AC) produces the same read as deleting 91-92 (CA).
    VCF spec requires the leftmost representation: anchor=G(pos90), deleted=CA -> REF=GCA, ALT=G.
    """
    row = _make_deletion_row(_FANC_AMPLICON, 92, 94)
    edits = list(vcf._edits_from_deletions(row, "FANC", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert chrom == "FANC"
    assert vcf_pos == 91, f"Expected left-normalized VCF pos 91, got {vcf_pos}"
    assert ref == "GCA", f"Expected REF=GCA, got {ref}"
    assert alt == "G", f"Expected ALT=G, got {alt}"


def test_deletion_left_normalization_one_base_shift():
    """Deletion (94,95) in FANC amplicon should left-normalize to VCF pos 93, REF=AC, ALT=A.

    Position 93 and 94 are both 'C', so deleting pos 94 is equivalent to deleting pos 93.
    VCF spec requires the leftmost: anchor=A(pos92), deleted=C(pos93) -> REF=AC, ALT=A.
    """
    row = _make_deletion_row(_FANC_AMPLICON, 94, 95)
    edits = list(vcf._edits_from_deletions(row, "FANC", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert chrom == "FANC"
    assert vcf_pos == 93, f"Expected left-normalized VCF pos 93, got {vcf_pos}"
    assert ref == "AC", f"Expected REF=AC, got {ref}"
    assert alt == "A", f"Expected ALT=A, got {alt}"


def test_deletion_no_normalization_needed():
    """Deletion (90,91) in FANC amplicon needs no left-normalization.

    Position 89=A, 90=G - distinct bases, so deletion at 90 is already leftmost.
    Expected VCF: anchor=A(pos89), deleted=G(pos90) -> POS=90, REF=AG, ALT=A.
    """
    row = _make_deletion_row(_FANC_AMPLICON, 90, 91)
    edits = list(vcf._edits_from_deletions(row, "FANC", 1))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert chrom == "FANC"
    assert vcf_pos == 90
    assert ref == "AG"
    assert alt == "A"


def test_deletion_at_start_with_pos_offset():
    """Deletion at start of amplicon with pos=100 should give VCF POS=100, not 99.

    Amplicon: GATTACA, delete first base G -> REF=GA, ALT=A.
    With pos=100, VCF POS should be 100 (position of the first deleted base).
    """
    row = _make_deletion_row("GATTACA", 0, 1)
    edits = list(vcf._edits_from_deletions(row, "chr1", 100))
    assert len(edits) == 1
    chrom, vcf_pos, ref, alt, reads = edits[0]
    assert chrom == "chr1"
    assert vcf_pos == 100, f"Expected VCF pos 100, got {vcf_pos}"
    assert ref == "GA", f"Expected REF=GA, got {ref}"
    assert alt == "A", f"Expected ALT=A, got {alt}"
