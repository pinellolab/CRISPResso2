from CRISPResso2 import CRISPRessoPooledCORE


def test_calculate_aligned_samtools_exclude_flags():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('0') == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_4():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('4') == hex(0x904)


def test_calculate_aligned_samtools_exclude_flags_9():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('9') == hex(0x909)


def test_calculate_aligned_samtools_exclude_flags_0x100():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('0x100') == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_010():
    """This tests for proper handling of octal numbers."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('010') == hex(0x908)


def test_calculate_aligned_samtools_exclude_flags_comma():
    """This tests for proper handling of commas."""
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('0,4') == f'0,4,{hex(0x900)}'
