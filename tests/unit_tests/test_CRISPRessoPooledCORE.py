from CRISPResso2 import CRISPRessoPooledCORE


def test_calculate_aligned_samtools_exclude_flags():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('0') == hex(0x900)


def test_calculate_aligned_samtools_exclude_flags_4():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('4') == hex(0x904)


def test_calculate_aligned_samtools_exclude_flags_9():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('9') == hex(0x909)


def test_calculate_aligned_samtools_exclude_flags_0x100():
    assert CRISPRessoPooledCORE.calculate_aligned_samtools_exclude_flags('0x100') == hex(0x900)
