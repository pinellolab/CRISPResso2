{
  "id": "dbd4df97",
  "title": "Fix multi-insertion bug: extract correct inserted bases using alignment index",
  "tags": [
    "bug",
    "vcf",
    "tdd"
  ],
  "status": "closed",
  "created_at": "2026-02-10T18:17:30.213Z"
}

## Goal

Fix `_edits_from_insertions` in `writers/vcf.py` to correctly extract inserted bases for reads with multiple insertions. Currently the second element of `insertion_coordinates` (a ref position) is used as an alignment index, which only works when there are no prior gaps in the alignment.

**Current behavior:** For reads with multiple insertions, the second and subsequent insertions extract wrong bases from the aligned sequence.
**Expected behavior:** All insertions extract the correct inserted bases regardless of prior gaps.

## Scope

- Fix inserted base extraction in `_edits_from_insertions`
- Rename misleading variables: `right_anchor_ref_pos` → `left_anchor_ref_pos`, `aligned_start` → `right_anchor_ref_pos`
- Keep the `== ref_len` defensive branch with corrected variable name and explanatory comment
- Update `design_docs/LEFT_NORMALIZATION.md` to correctly describe `insertion_coordinates` as `(left_anchor_ref_pos, right_anchor_ref_pos)`
- Add new tests only; don't modify existing tests

## Approach (TDD)

### Step 1: Write failing tests

**`test_build_alt_map.py`** — new test via `build_edit_counts`:
- Two insertions in the same read, verify both produce correct REF/ALT

**`test_write_vcf.py`** — new test via `_edits_from_insertions` directly:
- Row with two insertions, verify correct bases extracted for both

Use only DNA letters (A, C, T, G) in test sequences.

### Step 2: Verify tests fail (second insertion extracts wrong bases)

### Step 3: Fix `_edits_from_insertions`

- Rename variables to match actual semantics from `find_indels_substitutions`:
  - First element: `left_anchor_ref_pos` (was `right_anchor_ref_pos`)
  - Second element: `right_anchor_ref_pos` (was `aligned_start`)
- Replace base extraction:
  ```python
  # Old: ins_bases = aln_str[aligned_start:aligned_start + size]
  # New:
  right_anchor_idx = ref_positions.index(right_anchor_ref_pos)
  ins_bases = aln_str[right_anchor_idx - size:right_anchor_idx]
  ```
- Add comment on the defensive `left_anchor_ref_pos == ref_len` branch explaining it's currently unreachable because `find_indels_substitutions` doesn't record trailing insertions

### Step 4: Update `design_docs/LEFT_NORMALIZATION.md`

Fix the description of `insertion_coordinates` from `(anchor_ref_pos, aligned_start)` to `(left_anchor_ref_pos, right_anchor_ref_pos)` — both are reference positions.

### Step 5: Verify all tests pass (new + existing)

## Done when

- New multi-insertion tests pass
- Existing tests still pass
- Full `test_writers/` suite green
- Design doc accurately describes `insertion_coordinates`
