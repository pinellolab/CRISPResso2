"""Unit tests for CRISPRessoCompareCORE."""

import pytest

pytest.importorskip("scipy")

from CRISPResso2 import CRISPRessoCompareCORE

from copy import deepcopy
import pytest


@pytest.fixture(scope='function')
def run_info():
    return {
        'results': {
            'refs': {
                'Reference': {
                    'sequence':'CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG',
                    'sgRNA_orig_sequences': ['GGCCCTTAAAA'],
                    'sgRNA_cut_points': [50],
                    'allele_frequency_files': ['Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt'],
                },
            },
        },
    }


@pytest.fixture(scope='function')
def run_info_1(run_info):
    return deepcopy(run_info)


@pytest.fixture(scope='function')
def run_info_2(run_info):
    return deepcopy(run_info)


def test_get_matching_allele_files(run_info):
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info, run_info)
    assert matching_allele_files == [('Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt', 'Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt')]


def test_get_matching_allele_files_different_cut_points(run_info_1, run_info_2):
    run_info_2['results']['refs']['Reference']['sgRNA_cut_points'] = [50, 51]
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_different_guides(run_info_1, run_info_2):
    run_info_2['results']['refs']['Reference']['sgRNA_orig_sequences'] = ['GGCCCTTAAAC']
    run_info_2['results']['refs']['Reference']['allele_frequency_files'] = ['Alleles_frequency_table_around_sgRNA_GGCCCTTAAAC.txt']
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == []


def test_get_matching_allele_files_multiple_alleles(run_info_1, run_info_2):
    run_info_1['results']['refs']['Other_Amplicon'] = deepcopy(run_info_1['results']['refs']['Reference'])
    run_info_1['results']['refs']['Other_Amplicon']['sequence'] = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    run_info_1['results']['refs']['Other_Amplicon']['allele_frequency_files'] = ['Other_Amplicon.Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt']
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [('Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt', 'Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt')]


def test_get_matching_allele_files_different_amplicon_names_same_sequence(run_info_1, run_info_2):
    run_info_2['results']['refs']['Other_Amplicon'] = deepcopy(run_info_1['results']['refs']['Reference'])
    del run_info_2['results']['refs']['Reference']
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [('Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt', 'Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt')]


def test_get_matching_allele_files_some_different_guides(run_info_1, run_info_2):
    run_info_1['results']['refs']['Reference']['sgRNA_orig_sequences'] += ['AAAAAAAAAAAAAAAAAAA']
    run_info_1['results']['refs']['Reference']['allele_frequency_files'] += ['Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt']
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [('Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt', 'Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt')]


def test_get_matching_allele_files_multiple_guides(run_info_1, run_info_2):
    run_info_1['results']['refs']['Reference']['sgRNA_orig_sequences'] += ['AAAAAAAAAAAAAAAAAAA']
    run_info_1['results']['refs']['Reference']['allele_frequency_files'] += ['Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt']
    run_info_2['results']['refs']['Reference']['sgRNA_orig_sequences'] += ['AAAAAAAAAAAAAAAAAAA']
    run_info_2['results']['refs']['Reference']['allele_frequency_files'] += ['Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt']
    matching_allele_files = CRISPRessoCompareCORE.get_matching_allele_files(run_info_1, run_info_2)
    assert matching_allele_files == [
        ('Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt', 'Alleles_frequency_table_around_sgRNA_GGCCCTTAAAA.txt'),
        ('Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt', 'Alleles_frequency_table_around_sgRNA_AAAAAAAAAAAAAAAAAAA.txt'),
    ]
