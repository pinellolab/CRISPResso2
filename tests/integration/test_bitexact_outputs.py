import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DIR = REPO_ROOT / "tests"
EXPECTED_DIR = TESTS_DIR / "expectedResults"
HELPER = TESTS_DIR / "helpers" / "run_crispresso.py"

AMP_SEQ = (
    "CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAG"
    "CACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGT"
    "GGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG"
)
GUIDE_SEQ = "GGAATCCCTTCTGCAGCACC"
HDR_SEQ = (
    "CGGCCGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCTGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAG"
    "CTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTG"
    "CGCCACATCCATCGGCGCTTTGGTCGG"
)
CODING_SEQ = "GGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTT"


def _run_module(module_name, args, tmp_path, cwd, extra_env=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO_ROOT), env.get("PYTHONPATH", "")]
    ).strip(os.pathsep)
    env["MPLCONFIGDIR"] = str(tmp_path / "mplconfig")
    env["PYTHON"] = sys.executable
    scripts_dir = REPO_ROOT / "scripts"
    env["PATH"] = os.pathsep.join([str(scripts_dir), env.get("PATH", "")])
    if extra_env:
        env.update(extra_env)
    cmd = [sys.executable, str(HELPER), module_name, "--"] + args
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  cwd: {cwd}\n"
            f"  stdout:\n{result.stdout}\n"
            f"  stderr:\n{result.stderr}\n"
        )


def _assert_files_equal(actual_path, expected_path):
    actual_bytes = Path(actual_path).read_bytes()
    expected_bytes = Path(expected_path).read_bytes()
    assert actual_bytes == expected_bytes, f"Mismatch: {actual_path} != {expected_path}"


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_fanc_cas9_bitexact(tmp_path):
    out_dir = tmp_path / "fanc_cas9"
    args = [
        "-r1", str(TESTS_DIR / "FANC.Cas9.fastq"),
        "-a", AMP_SEQ,
        "-g", GUIDE_SEQ,
        "-p", "1",
        "--suppress_report",
        "--suppress_plots",
        "--output_folder", str(out_dir),
    ]
    _run_module("CRISPRessoCORE", args, tmp_path, TESTS_DIR)

    actual_root = out_dir / "CRISPResso_on_FANC.Cas9"
    expected_root = EXPECTED_DIR / "CRISPResso_on_FANC.Cas9"
    _assert_files_equal(
        actual_root / "Nucleotide_frequency_table.txt",
        expected_root / "Nucleotide_frequency_table.txt",
    )
    _assert_files_equal(
        actual_root / "CRISPResso_quantification_of_editing_frequency.txt",
        expected_root / "CRISPResso_quantification_of_editing_frequency.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_fanc_cas9_stream_align_bitexact(tmp_path):
    out_dir = tmp_path / "fanc_cas9_stream"
    args = [
        "-r1", str(TESTS_DIR / "FANC.Cas9.fastq"),
        "-a", AMP_SEQ,
        "-g", GUIDE_SEQ,
        "-p", "1",
        "--suppress_report",
        "--suppress_plots",
        "--output_folder", str(out_dir),
    ]
    extra_env = {"CRISPRESSO_FORCE_STREAM_ALIGN": "1"}
    _run_module("CRISPRessoCORE", args, tmp_path, TESTS_DIR, extra_env=extra_env)

    actual_root = out_dir / "CRISPResso_on_FANC.Cas9"
    expected_root = EXPECTED_DIR / "CRISPResso_on_FANC.Cas9"
    _assert_files_equal(
        actual_root / "Nucleotide_frequency_table.txt",
        expected_root / "Nucleotide_frequency_table.txt",
    )
    _assert_files_equal(
        actual_root / "CRISPResso_quantification_of_editing_frequency.txt",
        expected_root / "CRISPResso_quantification_of_editing_frequency.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_params_bitexact(tmp_path):
    out_dir = tmp_path / "params"
    args = [
        "-r1", str(TESTS_DIR / "FANC.Cas9.fastq"),
        "-a", AMP_SEQ,
        "-g", GUIDE_SEQ,
        "-e", HDR_SEQ,
        "-c", CODING_SEQ,
        "--dump",
        "-qwc", "20-30_45-50",
        "-q", "30",
        "--default_min_aln_score", "80",
        "-an", "FANC",
        "-n", "params",
        "--base_editor_output",
        "-fg", "AGCCTTGCAGTGGGCGCGCTA,CCCACTGAAGGCCC",
        "--dsODN", "GCTAGATTTCCCAAGAAGA",
        "-gn", "hi",
        "-fgn", "dear",
        "-p", "1",
        "--suppress_report",
        "--suppress_plots",
        "--output_folder", str(out_dir),
    ]
    _run_module("CRISPRessoCORE", args, tmp_path, TESTS_DIR)

    actual_root = out_dir / "CRISPResso_on_params"
    expected_root = EXPECTED_DIR / "CRISPResso_on_params"
    _assert_files_equal(
        actual_root / "FANC.Nucleotide_frequency_table.txt",
        expected_root / "FANC.Nucleotide_frequency_table.txt",
    )
    _assert_files_equal(
        actual_root / "CRISPResso_quantification_of_editing_frequency.txt",
        expected_root / "CRISPResso_quantification_of_editing_frequency.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_batch_bitexact(tmp_path):
    out_dir = tmp_path / "batch"
    args = [
        "-bs", str(TESTS_DIR / "FANC.local.batch"),
        "-a", AMP_SEQ,
        "-g", GUIDE_SEQ,
        "-p", "1",
        "--base_editor_output",
        "-n", "FANC",
        "--suppress_report",
        "--suppress_plots",
        "--batch_output_folder", str(out_dir),
    ]
    _run_module("CRISPRessoBatchCORE", args, tmp_path, TESTS_DIR)

    actual_root = out_dir / "CRISPRessoBatch_on_FANC"
    expected_root = EXPECTED_DIR / "CRISPRessoBatch_on_FANC"
    _assert_files_equal(
        actual_root / "MODIFICATION_FREQUENCY_SUMMARY.txt",
        expected_root / "MODIFICATION_FREQUENCY_SUMMARY.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_pooled_bitexact(tmp_path):
    out_dir = tmp_path / "pooled"
    args = [
        "-r1", str(TESTS_DIR / "Both.Cas9.fastq"),
        "-f", str(TESTS_DIR / "Cas9.amplicons.txt"),
        "-p", "1",
        "--keep_intermediate",
        "--min_reads_to_use_region", "100",
        "--suppress_report",
        "--suppress_plots",
        "--output_folder", str(out_dir),
    ]
    _run_module("CRISPRessoPooledCORE", args, tmp_path, TESTS_DIR)

    actual_root = out_dir / "CRISPRessoPooled_on_Both.Cas9"
    expected_root = EXPECTED_DIR / "CRISPRessoPooled_on_Both.Cas9"
    _assert_files_equal(
        actual_root / "SAMPLES_QUANTIFICATION_SUMMARY.txt",
        expected_root / "SAMPLES_QUANTIFICATION_SUMMARY.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_crispresso_wgs_bitexact(tmp_path):
    if shutil.which("samtools") is None:
        pytest.skip("samtools not available")

    out_dir = tmp_path / "wgs"
    args = [
        "-b", str(TESTS_DIR / "Both.Cas9.fastq.smallGenome.bam"),
        "-r", str(TESTS_DIR / "smallGenome" / "smallGenome.fa"),
        "-f", str(TESTS_DIR / "Cas9.regions.txt"),
        "-p", "1",
        "--suppress_report",
        "--suppress_plots",
        "--output_folder", str(out_dir),
    ]
    _run_module("CRISPRessoWGSCORE", args, tmp_path, TESTS_DIR)

    actual_root = out_dir / "CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome"
    expected_root = EXPECTED_DIR / "CRISPRessoWGS_on_Both.Cas9.fastq.smallGenome"
    _assert_files_equal(
        actual_root / "SAMPLES_QUANTIFICATION_SUMMARY.txt",
        expected_root / "SAMPLES_QUANTIFICATION_SUMMARY.txt",
    )


@pytest.mark.slow
@pytest.mark.integration
def test_parallel_parity_fanc_untreated(tmp_path):
    serial_dir = tmp_path / "serial"
    parallel_dir = tmp_path / "parallel"

    base_args = [
        "-r1", str(TESTS_DIR / "FANC.Untreated.fastq"),
        "-a", AMP_SEQ,
        "-g", GUIDE_SEQ,
        "--suppress_report",
        "--suppress_plots",
    ]

    _run_module(
        "CRISPRessoCORE",
        base_args + ["-p", "1", "--output_folder", str(serial_dir)],
        tmp_path,
        TESTS_DIR,
    )
    _run_module(
        "CRISPRessoCORE",
        base_args + ["-p", "2", "--output_folder", str(parallel_dir)],
        tmp_path,
        TESTS_DIR,
    )

    serial_root = serial_dir / "CRISPResso_on_FANC.Untreated"
    parallel_root = parallel_dir / "CRISPResso_on_FANC.Untreated"
    for filename in [
        "Nucleotide_frequency_table.txt",
        "CRISPResso_quantification_of_editing_frequency.txt",
    ]:
        _assert_files_equal(serial_root / filename, parallel_root / filename)
