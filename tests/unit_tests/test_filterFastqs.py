"""Tests for filterFastqs module."""

import gzip
import io
import os
import tempfile

import pytest

from CRISPResso2 import filterFastqs


def create_fastq_content(records):
    """Create FASTQ content string from records.

    Args:
        records: List of tuples (name, sequence, quality)

    Returns:
        String of FASTQ content
    """
    lines = []
    for name, seq, qual in records:
        lines.append(f"@{name}")
        lines.append(seq)
        lines.append("+")
        lines.append(qual)
    return "\n".join(lines) + "\n" if lines else ""


def write_fastq(filepath, records):
    """Write a FASTQ file from records."""
    with open(filepath, "w") as f:
        f.write(create_fastq_content(records))


def read_fastq_records(filepath):
    """Read FASTQ file and return list of (name, seq, qual) tuples."""
    records = []
    with open(filepath) as f:
        lines = f.readlines()
    for i in range(0, len(lines), 4):
        if i + 3 < len(lines):
            name = lines[i].strip().lstrip("@")
            seq = lines[i + 1].strip()
            qual = lines[i + 3].strip()
            records.append((name, seq, qual))
    return records


# =============================================================================
# Tests for run_mBPN (convert low quality bases to N)
# =============================================================================


def test_run_mBPN_high_quality_unchanged():
    """Test that high quality bases are not converted to N."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Quality 'I' = 40 (Phred+33), well above threshold
        write_fastq(in_file, [("read1", "ATCG", "IIII")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBPN(f_in, f_out, None, None, 20)

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][1] == "ATCG"  # Sequence unchanged


def test_run_mBPN_low_quality_converted():
    """Test that low quality bases are converted to N."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Quality '!' = 0 (Phred+33), below threshold of 20
        write_fastq(in_file, [("read1", "ATCG", "!!!!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBPN(f_in, f_out, None, None, 20)

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][1] == "NNNN"  # All converted to N


def test_run_mBPN_mixed_quality():
    """Test that only low quality bases are converted to N."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # 'I' = 40, '!' = 0, threshold = 20
        # Positions 0 and 2 should stay, 1 and 3 should become N
        write_fastq(in_file, [("read1", "ATCG", "I!I!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBPN(f_in, f_out, None, None, 20)

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][1] == "ANCN"


def test_run_mBPN_multiple_reads():
    """Test run_mBPN with multiple reads."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [
            ("read1", "ATCG", "IIII"),  # High quality
            ("read2", "GCTA", "!!!!"),  # Low quality
        ])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBPN(f_in, f_out, None, None, 20)

        records = read_fastq_records(out_file)
        assert len(records) == 2
        assert records[0][1] == "ATCG"
        assert records[1][1] == "NNNN"


def test_run_mBPN_empty_file():
    """Test run_mBPN with empty input file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBPN(f_in, f_out, None, None, 20)

        assert os.path.getsize(out_file) == 0


# =============================================================================
# Tests for run_mRQ (filter by average read quality)
# =============================================================================


def test_run_mRQ_high_quality_passes():
    """Test that reads with high average quality pass."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Quality 'I' = 40, average = 40, threshold = 20
        write_fastq(in_file, [("read1", "ATCG", "IIII")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ(f_in, f_out, None, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 1


def test_run_mRQ_low_quality_filtered():
    """Test that reads with low average quality are filtered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Quality '!' = 0, average = 0, threshold = 20
        write_fastq(in_file, [("read1", "ATCG", "!!!!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ(f_in, f_out, None, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 0


def test_run_mRQ_borderline_quality():
    """Test that reads at exactly the threshold pass."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Quality '5' = 20 (Phred+33), average = 20, threshold = 20
        write_fastq(in_file, [("read1", "ATCG", "5555")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ(f_in, f_out, None, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 1


def test_run_mRQ_mixed_reads():
    """Test run_mRQ with mixed quality reads."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [
            ("read1", "ATCG", "IIII"),  # High quality, should pass
            ("read2", "GCTA", "!!!!"),  # Low quality, should fail
            ("read3", "AAAA", "IIII"),  # High quality, should pass
        ])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ(f_in, f_out, None, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 2
        assert records[0][0] == "read1"
        assert records[1][0] == "read3"


# =============================================================================
# Tests for run_mBP (filter by minimum base pair quality)
# =============================================================================


def test_run_mBP_all_high_quality_passes():
    """Test that reads with all high quality bases pass."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # All bases quality 40, threshold 20
        write_fastq(in_file, [("read1", "ATCG", "IIII")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP(f_in, f_out, 20, None, None)

        records = read_fastq_records(out_file)
        assert len(records) == 1


def test_run_mBP_any_low_quality_filtered():
    """Test that reads with any low quality base are filtered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # One low quality base at position 1
        write_fastq(in_file, [("read1", "ATCG", "I!II")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP(f_in, f_out, 20, None, None)

        records = read_fastq_records(out_file)
        assert len(records) == 0


def test_run_mBP_all_low_quality_filtered():
    """Test that reads with all low quality bases are filtered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [("read1", "ATCG", "!!!!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP(f_in, f_out, 20, None, None)

        records = read_fastq_records(out_file)
        assert len(records) == 0


# =============================================================================
# Tests for run_mBP_mRQ (combined filters)
# =============================================================================


def test_run_mBP_mRQ_both_pass():
    """Test that reads passing both filters are kept."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [("read1", "ATCG", "IIII")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP_mRQ(f_in, f_out, 20, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 1


def test_run_mBP_mRQ_fail_mBP():
    """Test that reads failing mBP are filtered even if mRQ passes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # One low quality base, but average still passes
        # 'I' = 40, '!' = 0, average = 30
        write_fastq(in_file, [("read1", "ATCG", "III!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP_mRQ(f_in, f_out, 20, 20, None)

        records = read_fastq_records(out_file)
        assert len(records) == 0


# =============================================================================
# Tests for mBP filter via filterFastqs main function
# =============================================================================


def test_filterFastqs_mBP_filters_low_quality():
    """Test that reads with any base below mBP threshold are filtered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # One base below mBP threshold
        write_fastq(in_file, [("read1", "ATCG", "!III")])

        filterFastqs.filterFastqs(
            fastq_r1=in_file,
            fastq_r1_out=out_file,
            min_bp_qual_in_read=20,
        )

        records = read_fastq_records(out_file)
        assert len(records) == 0


# =============================================================================
# Tests for run_mRQ_mBPN (combined filters)
# =============================================================================


def test_run_mRQ_mBPN_passes_and_converts():
    """Test mRQ filtering and mBPN conversion together."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Average quality passes mRQ, some bases below mBPN
        # 'I' = 40, '5' = 20, average = 30, mRQ threshold = 25
        write_fastq(in_file, [("read1", "ATCG", "II55")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ_mBPN(f_in, f_out, None, 25, 30)

        records = read_fastq_records(out_file)
        assert len(records) == 1
        # Last two bases (qual 20) should be N, first two (qual 40) kept
        assert records[0][1] == "ATNN"


def test_run_mRQ_mBPN_fails_mRQ():
    """Test that reads failing mRQ are filtered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Average quality fails
        write_fastq(in_file, [("read1", "ATCG", "!!!!")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mRQ_mBPN(f_in, f_out, None, 20, 30)

        records = read_fastq_records(out_file)
        assert len(records) == 0


# =============================================================================
# Tests for run_mBP_mRQ_mBPN (all three filters combined)
# =============================================================================


def test_run_mBP_mRQ_mBPN_all_pass():
    """Test that reads passing all filters are kept and converted."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        # All bases above mBP (10), average above mRQ (20), some below mBPN (35)
        # '5' = 20, 'I' = 40
        write_fastq(in_file, [("read1", "ATCG", "55II")])

        with open(in_file, "rb") as f_in, open(out_file, "w") as f_out:
            filterFastqs.run_mBP_mRQ_mBPN(f_in, f_out, 10, 20, 35)

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][1] == "NNCG"


# =============================================================================
# Tests for filterFastqs main function
# =============================================================================


def test_filterFastqs_single_file_mBPN():
    """Test filterFastqs with single file and mBPN filter."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [
            ("read1", "ATCG", "IIII"),
            ("read2", "GCTA", "!!!!"),
        ])

        filterFastqs.filterFastqs(
            fastq_r1=in_file,
            fastq_r1_out=out_file,
            min_bp_qual_or_N=20,
        )

        records = read_fastq_records(out_file)
        assert len(records) == 2
        assert records[0][1] == "ATCG"
        assert records[1][1] == "NNNN"


def test_filterFastqs_single_file_mRQ():
    """Test filterFastqs with single file and mRQ filter."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [
            ("read1", "ATCG", "IIII"),
            ("read2", "GCTA", "!!!!"),
        ])

        filterFastqs.filterFastqs(
            fastq_r1=in_file,
            fastq_r1_out=out_file,
            min_av_read_qual=20,
        )

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][0] == "read1"


def test_filterFastqs_single_file_mBP():
    """Test filterFastqs with single file and mBP filter."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq")
        out_file = os.path.join(tmpdir, "output.fastq")

        write_fastq(in_file, [
            ("read1", "ATCG", "IIII"),
            ("read2", "GCTA", "I!II"),  # One low quality base
        ])

        filterFastqs.filterFastqs(
            fastq_r1=in_file,
            fastq_r1_out=out_file,
            min_bp_qual_in_read=20,
        )

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][0] == "read1"


def test_filterFastqs_nonexistent_file():
    """Test filterFastqs raises exception for nonexistent file."""
    with pytest.raises(Exception, match="does not exist"):
        filterFastqs.filterFastqs(
            fastq_r1="/nonexistent/path/file.fastq",
            min_bp_qual_or_N=20,
        )


def test_filterFastqs_gzipped_input():
    """Test filterFastqs with gzipped input file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "input.fastq.gz")
        out_file = os.path.join(tmpdir, "output.fastq")

        # Create gzipped FASTQ
        content = create_fastq_content([
            ("read1", "ATCG", "IIII"),
            ("read2", "GCTA", "!!!!"),
        ])
        with gzip.open(in_file, "wt") as f:
            f.write(content)

        filterFastqs.filterFastqs(
            fastq_r1=in_file,
            fastq_r1_out=out_file,
            min_av_read_qual=20,
        )

        records = read_fastq_records(out_file)
        assert len(records) == 1
        assert records[0][0] == "read1"
