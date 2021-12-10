from CRISPResso2 import CRISPRessoCOREResources

def test_find_indels_substitutions():
    # no insertion or deletion
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        [18, 19],
    )
    assert payload['insertion_n'] == 0
    assert payload['deletion_n'] == 0
    assert payload['substitution_n'] == 0

    # deletion outside of quantification window
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG',
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        [21, 22],
    )
    assert payload['insertion_n'] == 0
    assert payload['substitution_n'] == 0
    assert payload['deletion_n'] == 0

    # deletion overlap quantification window
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG',
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        [18, 19],
    )
    assert payload['insertion_n'] == 0
    assert payload['substitution_n'] == 0
    assert payload['deletion_n'] == 3
    assert payload['all_deletion_positions'] == [18, 19, 20]
    assert payload['deletion_positions'] == [18, 19, 20]
    assert payload['deletion_coordinates'] == [(18, 21)]

    # insertion outside of quantification window
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        'CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG',
        [21, 22],
    )
    assert payload['insertion_n'] == 0
    assert payload['substitution_n'] == 0
    assert payload['deletion_n'] == 0

    # insertion overlap quantification window
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        'CATGGAATCCCTTCTGCA---CCTGGATCGCTTTTCCGAG',
        [18, 19],
    )
    assert payload['insertion_n'] == 0
    assert payload['all_insertion_positions'] == [17, 18]
    assert payload['all_insertion_left_positions'] == [17]
    assert payload['deletion_n'] == 0
    assert payload['substitution_n'] == 0

    # deletion at ends
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        '---GGAATCCCTTCTGCAGCACCTGGATCGCTTTTC----',
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        [18, 19],
    )
    assert payload['insertion_n'] == 0
    assert payload['deletion_n'] == 0
    assert payload['substitution_n'] == 0

    # insertion at ends
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'CATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAG',
        '----GAATCCCTTCTGCAGCACCTGGATCGCTTTTC----',
        [18, 19],
    )
    assert payload['insertion_n'] == 0
    assert payload['deletion_n'] == 0
    assert payload['substitution_n'] == 0

    # insertion at the edge of include_indx_set
    payload = CRISPRessoCOREResources.find_indels_substitutions(
        'TAATCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAAGAGGGCGGCTTTGGGCGGGGTC-CAGTTCCGGGATTA--GCGAACTTAGAGCAC-----ACGTCTGAACTCCAGTCACCGATGTATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACTTACTCTCACTTAACTCTTGCTTCCCTCCTGACGCCGATG',
        '----CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGC------------ACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCA-AGCACTACCTACGTCAGCACCTGGGACCCCGCCAC------CGTGCGCCGGGC----CTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGG--CGCTTTGGTCGG-----',
        [91, 92],
    )
    correct_payload = {
        'all_insertion_positions': [91, 92, 126, 127, 161, 162, 173, 174, 210, 211],
        'all_insertion_left_positions': [91, 126, 161, 173, 210],
        'insertion_positions': [91, 92],
        'insertion_coordinates': [(91, 92)],
        'insertion_sizes': [12],
        'insertion_n': 12,
        'all_deletion_positions': [101, 116, 117, 132, 133, 134, 135, 136],
        'deletion_positions': [],
        'deletion_coordinates': [],
        'deletion_sizes': [],
        'deletion_n': 0.0,
        'all_substitution_positions': [90, 91, 92, 93, 95, 98, 103, 104, 110, 112, 115, 121, 122, 125, 142, 144, 147, 148, 149, 150, 152, 154, 158, 159, 160, 161, 165, 166, 171, 172, 173, 178, 180, 181, 182, 183, 184, 185, 186, 187, 188, 191, 192, 195, 196, 197, 198, 200, 201, 203, 205, 206, 208, 210, 212, 215, 216, 217, 219, 222],
        'substitution_positions': [91, 92],
        'substitution_n': 2,
        'ref_positions': [-1, -1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, -92, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, -127, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, -162, -162, -162, -162, -162, -162, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, -174, -174, -174, -174, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, -211, -211, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, -223, -223, -223, -223, -223],
    }
    assert payload['all_insertion_positions'] == correct_payload['all_insertion_positions']
    assert payload['all_insertion_left_positions'] == correct_payload['all_insertion_left_positions']
    assert payload['insertion_positions'] == correct_payload['insertion_positions']
    assert payload['insertion_coordinates'] == correct_payload['insertion_coordinates']
    assert payload['insertion_sizes'] == correct_payload['insertion_sizes']
    assert payload['insertion_n'] == correct_payload['insertion_n']
    assert payload['all_deletion_positions'] == correct_payload['all_deletion_positions']
    assert payload['deletion_positions'] == correct_payload['deletion_positions']
    assert payload['deletion_sizes'] == correct_payload['deletion_sizes']
    assert payload['deletion_n'] == correct_payload['deletion_n']
    assert payload['all_substitution_positions'] == correct_payload['all_substitution_positions']
    assert payload['substitution_positions'] == correct_payload['substitution_positions']
    assert payload['substitution_n'] == correct_payload['substitution_n']
    assert payload['ref_positions'] == correct_payload['ref_positions']


if __name__ == "__main__":
    # execute only if run as a script
    test_find_indels_substitutions()
