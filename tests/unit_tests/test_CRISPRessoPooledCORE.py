from CRISPResso2 import CRISPRessoPooledCORE


def test_calculate_aligned_samtools_exclude_flags():
    """Test calculate_aligned_samtools_exclude_flags with flag 0."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("0") == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_4():
    """Test calculate_aligned_samtools_exclude_flags with flag 4."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("4") == hex(0x904)


def test_calculate_aligned_samtools_exclude_flags_9():
    """Test calculate_aligned_samtools_exclude_flags with flag 9."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("9") == hex(0x909)


def test_calculate_aligned_samtools_exclude_flags_0x100():
    """Test calculate_aligned_samtools_exclude_flags with flag 0x100."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("0x100") == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_010():
    """This tests for proper handling of octal numbers."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("010") == hex(0x908)


def test_calculate_aligned_samtools_exclude_flags_comma():
    """This tests for proper handling of commas."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("0,4") == f"0,4,{hex(0x900)}"


# =============================================================================
# Additional edge case tests
# =============================================================================


def test_calculate_aligned_samtools_exclude_flags_large():
    """Test with large flag value."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("4096")
    assert result == "0x1900"


def test_calculate_aligned_samtools_exclude_flags_multiple_commas():
    """Test with multiple comma-separated values."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("0,4,8")
    assert result == "0,4,8,0x900"


def test_calculate_aligned_samtools_exclude_flags_hex_lowercase():
    """Test with lowercase hex prefix."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("0x10")
    assert result == "0x910"


def test_calculate_aligned_samtools_exclude_flags_1():
    """Test with flag 1."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("1")
    assert result == hex(0x901)


def test_calculate_aligned_samtools_exclude_flags_2():
    """Test with flag 2."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("2")
    assert result == hex(0x902)


def test_calculate_aligned_samtools_exclude_flags_8():
    """Test with flag 8."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("8")
    assert result == hex(0x908)


def test_calculate_aligned_samtools_exclude_flags_16():
    """Test with flag 16."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("16")
    assert result == hex(0x910)


def test_calculate_aligned_samtools_exclude_flags_256():
    """Test with flag 256 (0x100) - secondary alignment."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("256")
    # 256 is already in base flags, so result should be same as 0
    assert result == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_2048():
    """Test with flag 2048 (0x800) - supplementary alignment."""
    result = CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags("2048")
    # 2048 is already in base flags, so result should be same as 0
    assert result == hex(0x900)
