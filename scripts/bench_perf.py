#!/usr/bin/env python3
import argparse
import os
import sys
import time
from types import SimpleNamespace

import importlib.metadata as _importlib_metadata
import numpy as np

# Ensure we import CRISPResso2 from the repo root when run from subdirs.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Avoid importlib.metadata failure when CRISPResso2 isn't installed as a package.
_orig_version = _importlib_metadata.version


def _version_patch(name):
    if name == "CRISPResso2":
        return "dev"
    return _orig_version(name)


_importlib_metadata.version = _version_patch

from CRISPResso2 import CRISPResso2Align
from CRISPResso2 import CRISPRessoCOREResources
from CRISPResso2 import CRISPRessoShared
from CRISPResso2 import CRISPRessoCORE


ALPHABET = np.array(list("ACGTN"), dtype="U1")


def bench(name, fn, iters):
    t0 = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t0
    rate = iters / dt if dt else float("inf")
    print(f"{name}\t{dt:.4f}s\t{rate:.1f} ops/s")
    return out


def make_random_seqs(rng, count, length):
    arr = rng.choice(ALPHABET, size=(count, length))
    return ["".join(row) for row in arr]


def make_aligned_pairs(rng, count, ref_len, ins_rate=0.02, del_rate=0.02, sub_rate=0.02):
    pairs = []
    ref = "".join(rng.choice(ALPHABET, size=ref_len))
    for _ in range(count):
        ref_aln = []
        read_aln = []
        for base in ref:
            if rng.random() < ins_rate:
                ref_aln.append("-")
                read_aln.append(rng.choice(ALPHABET))
            if rng.random() < del_rate:
                ref_aln.append(base)
                read_aln.append("-")
            else:
                if rng.random() < sub_rate:
                    read_base = rng.choice(ALPHABET)
                else:
                    read_base = base
                ref_aln.append(base)
                read_aln.append(read_base)
        pairs.append(("".join(read_aln), "".join(ref_aln)))
    return pairs


def build_refs(ref_seqs):
    refs = {}
    for i, seq in enumerate(ref_seqs):
        ref_name = f"ref{i+1}"
        refs[ref_name] = {
            "sequence": seq,
            "sequence_length": len(seq),
            "min_aln_score": -1,
            "gap_incentive": np.zeros(len(seq) + 1, dtype=int),
            "fw_seeds": [],
            "rc_seeds": [],
            "include_idxs": np.arange(len(seq)),
        }
    return refs


def build_args():
    return SimpleNamespace(
        aln_seed_count=0,
        aln_seed_min=0,
        needleman_wunsch_gap_open=-1,
        needleman_wunsch_gap_extend=-1,
        needleman_wunsch_bandwidth=-1,
        use_legacy_insertion_quantification=False,
        ignore_deletions=False,
        ignore_insertions=False,
        ignore_substitutions=False,
        assign_ambiguous_alignments_to_first_reference=False,
        expand_ambiguous_alignments=False,
        prime_editing_pegRNA_scaffold_seq="",
    )


def main():
    parser = argparse.ArgumentParser(description="CRISPResso2 micro-benchmarks")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--rc-seqs", type=int, default=20000)
    parser.add_argument("--rc-len", type=int, default=150)
    parser.add_argument("--indel-pairs", type=int, default=2000)
    parser.add_argument("--indel-len", type=int, default=200)
    parser.add_argument("--align-reads", type=int, default=200)
    parser.add_argument("--align-len", type=int, default=150)
    parser.add_argument("--align-refs", type=int, default=2)
    parser.add_argument("--no-align", action="store_true")
    args = parser.parse_args()

    print(f"CRISPResso2 version: {CRISPRessoShared.__version__}")
    print(f"Python: {sys.version.split()[0]}  NumPy: {np.__version__}")

    rng = np.random.default_rng(args.seed)

    rc_seqs = make_random_seqs(rng, args.rc_seqs, args.rc_len)
    bench(
        "reverse_complement",
        lambda: [CRISPRessoShared.reverse_complement(s) for s in rc_seqs],
        args.rc_seqs,
    )

    include_idxs = np.arange(args.indel_len)
    aligned_pairs = make_aligned_pairs(rng, max(1, args.indel_pairs), args.indel_len)
    bench(
        "find_indels_substitutions",
        lambda: [
            CRISPRessoCOREResources.find_indels_substitutions(r, a, include_idxs)
            for r, a in aligned_pairs
        ],
        len(aligned_pairs),
    )

    if not args.no_align:
        ref_seqs = make_random_seqs(rng, args.align_refs, args.align_len)
        refs = build_refs(ref_seqs)
        ref_names = list(refs.keys())
        aln_matrix = CRISPResso2Align.make_matrix()
        bench_reads = make_random_seqs(rng, args.align_reads, args.align_len)
        bench_args = build_args()
        bench(
            "get_new_variant_object",
            lambda: [
                CRISPRessoCORE.get_new_variant_object(
                    bench_args,
                    s,
                    refs,
                    ref_names,
                    aln_matrix,
                    (0, None),
                )
                for s in bench_reads
            ],
            len(bench_reads),
        )

        bench_args_banded = build_args()
        bench_args_banded.needleman_wunsch_bandwidth = 20
        bench(
            "get_new_variant_object_banded",
            lambda: [
                CRISPRessoCORE.get_new_variant_object(
                    bench_args_banded,
                    s,
                    refs,
                    ref_names,
                    aln_matrix,
                    (0, None),
                )
                for s in bench_reads
            ],
            len(bench_reads),
        )


if __name__ == "__main__":
    main()
