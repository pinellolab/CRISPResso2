"""Unit tests for CRISPRessoPlot.py"""

import pytest
from CRISPResso2 import CRISPRessoPlot

def test_remove_deletions_for_amino_acids():

    # no complete codons removed
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('A--T-C--G-T') == 'ATCGT'
    # one complete codon removed
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG') == 'CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG'
    # frameshift
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('CATGGAATCCCTTCTGC---CCTGGATCGCTTTTCCGAG') == 'CATGGAATCCCTTCTGCCCTGGATCGCTTTTCCGAG'
    # frameshift and complete codon removed
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('C--GGAATCCCTTCTGC---CCTGGATCGCTTTTCCGAG') == 'CGGAATCCCTTCTGC---CCTGGATCGCTTTTCCGAG'

    assert CRISPRessoPlot.remove_deletions_for_amino_acids('--CCC---GT----ACAT--') == 'CCC---GTACAT'
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('--CCC--GTC----ACAT--') == 'CCCGTC---ACAT'
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('--------------------') == '------------------'
    assert CRISPRessoPlot.remove_deletions_for_amino_acids('CATGGAATCCCTTCTG----ACCTGGATCGCTTTTCCGAG') == 'CATGGAATCCCTTCTGACCTGGATCGCTTTTCCGAG'

