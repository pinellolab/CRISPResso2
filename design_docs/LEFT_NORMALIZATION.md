# VCF Left-Normalization in `writers/vcf.py`

## Problem

CRISPResso2's aligner may **right-shift indels** in repetitive regions of the
reference sequence. For example, in `...GCACCT...` (positions 90-95 of the FANC
amplicon), deleting the `CA` at positions 91-92 produces the same read as
deleting `AC` at positions 92-93. The aligner chose 92-93, but the VCF spec
requires the **leftmost** equivalent representation (91-92).

The same applies to insertions: inserting a `C` into a `CCC` homopolymer run
can be placed at any position within the run. VCF requires the leftmost.

## Root Cause

`_edits_from_deletions` and `_edits_from_insertions` in `writers/vcf.py` were
passing aligner-provided `deletion_coordinates` and `insertion_coordinates`
directly into VCF records without left-normalizing them. The coordinates come
from `CRISPRessoCOREResources.find_indels_substitutions`, which assigns
positions based on alignment gaps without enforcing leftmost placement.

## How It Was Found

The `vcf-basic` integration test (in the CRISPResso2_tests repo worktree at
`cli_integration_tests/`) produced a VCF with two wrong deletion records:

| Expected | Actual (buggy) | Cause |
|----------|----------------|-------|
| `FANC 91 GCA G` | `FANC 92 CAC C` | 2bp deletion right-shifted by 1 |
| `FANC 93 AC A` | `FANC 94 CC C` | 1bp deletion right-shifted by 1 |

Both cases involve repetitive `C` bases where the aligner placed the deletion
one position to the right of the canonical leftmost position.

## Fix

Two functions were added to `writers/vcf.py`:

### `_left_normalize_deletion(start, end, ref_positions, ref_str)`

Shifts the half-open deletion window `[start, end)` left while the base just
before the deletion equals the last deleted base:

```python
while start > 0:
    if ref_str[ref_positions.index(start - 1)] == ref_str[ref_positions.index(end - 1)]:
        start -= 1
        end -= 1
    else:
        break
```

### `_left_normalize_insertion(anchor_ref_pos, ins_bases, ref_positions, ref_str)`

Rotates the inserted bases and shifts the anchor position left while the last
inserted base equals the base at the current anchor:

```python
while anchor_ref_pos > 0:
    if ins_bases[-1] == ref_str[ref_positions.index(anchor_ref_pos)]:
        ins_bases = ref_str[ref_positions.index(anchor_ref_pos)] + ins_bases[:-1]
        anchor_ref_pos -= 1
    else:
        break
```

Both functions use `ref_positions` (the alignment-index-to-reference-position
mapping) and `ref_str` (the aligned reference sequence) to look up bases,
making them correct even when the alignment contains other indels.

## Key Data Structures

- **`deletion_coordinates`**: List of `(start, end)` tuples, 0-based half-open
  `[start, end)`. Produced by `CRISPRessoCOREResources.find_indels_substitutions`.
- **`insertion_coordinates`**: List of `(left_anchor_ref_pos, right_anchor_ref_pos)`
  tuples. Both are 0-based reference positions produced by
  `CRISPRessoCOREResources.find_indels_substitutions`. `left_anchor_ref_pos` is
  the base before the insertion; `right_anchor_ref_pos` is the base after.
  To find the inserted bases in the aligned sequence, convert
  `right_anchor_ref_pos` to an alignment index via `ref_positions.index()` and
  take the `size` characters immediately before it.
- **`ref_positions`**: List mapping each alignment column index to its reference
  position (or -1 for inserted bases).

## Tests

Unit tests in `tests/unit_tests/test_writers/test_vcf/test_write_vcf.py`:

- `TestLeftNormalizeDeletion` (9 tests): direct tests of the deletion normalizer
- `TestLeftNormalizeInsertion` (8 tests): direct tests of the insertion normalizer
- `test_deletion_left_normalization_*` (3 tests): end-to-end through `_edits_from_deletions`
- `test_insertion_left_normalization_*` (3 tests): end-to-end through `_edits_from_insertions`

Integration test: `make vcf-basic test` in the CRISPResso2_tests worktree.
