"""Shared pytest fixtures for CRISPResso2 unit tests."""

import os
import tempfile

import pytest


@pytest.fixture
def temp_dir():
    """Provide a temporary directory that's cleaned up after tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_fastq(temp_dir):
    """Create a sample FASTQ file for testing."""
    filepath = os.path.join(temp_dir, "sample.fastq")
    with open(filepath, "w") as f:
        f.write("@read1\nATCGATCG\n+\nIIIIIIII\n")
        f.write("@read2\nGCTAGCTA\n+\nIIIIIIII\n")
    return filepath


@pytest.fixture
def sample_fastq_low_quality(temp_dir):
    """Create a sample FASTQ file with low quality scores."""
    filepath = os.path.join(temp_dir, "low_quality.fastq")
    with open(filepath, "w") as f:
        f.write("@read1\nATCGATCG\n+\n!!!!!!!!\n")  # Quality 0
        f.write("@read2\nGCTAGCTA\n+\n########\n")  # Quality 2
    return filepath


@pytest.fixture
def sample_fastq_mixed_quality(temp_dir):
    """Create a sample FASTQ file with mixed quality scores."""
    filepath = os.path.join(temp_dir, "mixed_quality.fastq")
    with open(filepath, "w") as f:
        f.write("@read1\nATCGATCG\n+\nIIIIIIII\n")  # High quality
        f.write("@read2\nGCTAGCTA\n+\n!!!!!!!!\n")  # Low quality
        f.write("@read3\nAAAAAAAA\n+\nIIIIIIII\n")  # High quality
    return filepath


@pytest.fixture
def empty_fastq(temp_dir):
    """Create an empty FASTQ file."""
    filepath = os.path.join(temp_dir, "empty.fastq")
    with open(filepath, "w") as f:
        f.close()
    return filepath


@pytest.fixture
def aln_matrix():
    """Load the EDNAFULL alignment matrix."""
    from CRISPResso2 import CRISPResso2Align

    return CRISPResso2Align.read_matrix("./CRISPResso2/EDNAFULL")


@pytest.fixture
def blosum62_matrix():
    """Load the BLOSUM62 alignment matrix."""
    from CRISPResso2 import CRISPResso2Align

    return CRISPResso2Align.read_matrix("./CRISPResso2/BLOSUM62")


def create_test_fastq(filepath, records):
    """Helper function to create test FASTQ files.

    Args:
        filepath: Path to create the file at
        records: List of tuples (name, sequence, quality)
    """
    with open(filepath, "w") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")
    return filepath
