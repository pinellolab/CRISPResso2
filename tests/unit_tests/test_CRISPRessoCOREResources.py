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


if __name__ == "__main__":
# execute only if run as a script
    test_find_indels_substitutions()
