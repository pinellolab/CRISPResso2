"""Tests for CRISPRessoMultiProcessing module."""

import os

import pandas as pd
import pytest

from CRISPResso2 import CRISPRessoMultiProcessing


# Module-level helper functions for multi-process tests.
# multiprocessing.Pool requires picklable (i.e. top-level) functions.
def _mp_square_all(arr):
    return [x ** 2 for x in arr]


def _mp_uppercase_all(arr):
    return [s.upper() for s in arr]


def _mp_identity(arr):
    return list(arr)


def _mp_double_column(df_chunk):
    df_chunk = df_chunk.copy()
    df_chunk["doubled"] = df_chunk["value"] * 2
    return df_chunk


def _mp_add_column(df_chunk):
    df_chunk = df_chunk.copy()
    df_chunk["new_col"] = "added"
    return df_chunk


# =============================================================================
# Tests for get_max_processes
# =============================================================================


def test_get_max_processes():
    """Test that get_max_processes returns os.cpu_count()."""
    max_procs = CRISPRessoMultiProcessing.get_max_processes()
    assert max_procs == os.cpu_count()


# =============================================================================
# Tests for wrapper function
# =============================================================================


def test_wrapper_returns_index_and_result():
    """Test that wrapper returns tuple of (index, func(arg))."""
    def add_one(x):
        return x + 1

    result = CRISPRessoMultiProcessing.wrapper(add_one, (0, 5))
    assert result == (0, 6)


def test_wrapper_with_none_result():
    """Test wrapper handles None result from function."""
    def return_none(x):
        return None

    result = CRISPRessoMultiProcessing.wrapper(return_none, (0, "anything"))
    assert result == (0, None)


def test_wrapper_with_complex_return():
    """Test wrapper handles complex return types."""
    def return_dict(x):
        return {"input": x, "processed": True}

    result = CRISPRessoMultiProcessing.wrapper(return_dict, (5, "data"))
    assert result == (5, {"input": "data", "processed": True})


# =============================================================================
# Tests for run_subprocess
# =============================================================================


def test_run_subprocess_success():
    """Test run_subprocess with successful command returns 0."""
    assert CRISPRessoMultiProcessing.run_subprocess("true") == 0


def test_run_subprocess_failure():
    """Test run_subprocess with failing command returns non-zero."""
    assert CRISPRessoMultiProcessing.run_subprocess("false") != 0


def test_run_subprocess_exit_code():
    """Test run_subprocess captures specific exit codes."""
    assert CRISPRessoMultiProcessing.run_subprocess("sh -c 'exit 5'") == 5


# =============================================================================
# Tests for run_function_on_array_chunk_parallel — single process
# =============================================================================


def test_run_function_on_array_chunk_parallel_single_process():
    """Test single-process path calls function directly on full array."""
    def square_all(arr):
        return [x ** 2 for x in arr]

    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        [1, 2, 3, 4, 5], square_all, n_processes=1
    )
    assert result == [1, 4, 9, 16, 25]


def test_run_function_on_array_chunk_parallel_empty_array():
    """Test single-process path with empty array."""
    def identity(arr):
        return arr

    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        [], identity, n_processes=1
    )
    assert result == []


def test_run_function_on_array_chunk_parallel_exception_handling():
    """Test that exceptions in the function are propagated."""
    def raise_error(arr):
        raise ValueError("Test error")

    with pytest.raises(ValueError, match="Test error"):
        CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
            [1, 2, 3], raise_error, n_processes=1
        )


# =============================================================================
# Tests for run_function_on_array_chunk_parallel — multi process
# =============================================================================


def test_run_function_on_array_chunk_parallel_multi_process():
    """Test multi-process path splits array into chunks and flattens results.

    Uses 25 elements with n_processes=2 so chunk size is max(10, 12)=12,
    producing 2 real chunks that exercise pool.map_async and result flattening.
    """
    input_array = list(range(25))
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, _mp_square_all, n_processes=2
    )
    assert result == [x ** 2 for x in range(25)]


def test_run_function_on_array_chunk_parallel_multi_process_strings():
    """Test multi-process path with string data (pickling of non-numeric types).

    Uses 20 elements so chunk size max(10, 10)=10 produces 2 chunks.
    """
    input_array = [f"word_{i}" for i in range(20)]
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        input_array, _mp_uppercase_all, n_processes=2
    )
    assert result == [f"WORD_{i}" for i in range(20)]


def test_run_function_on_array_chunk_parallel_multi_process_empty():
    """Test multi-process path with empty array doesn't error."""
    result = CRISPRessoMultiProcessing.run_function_on_array_chunk_parallel(
        [], _mp_identity, n_processes=2
    )
    assert result == []


# =============================================================================
# Tests for run_pandas_apply_parallel — single process
# =============================================================================


def test_run_pandas_apply_parallel_single_process():
    """Test single-process path calls function directly on full dataframe."""
    def double_column(df_chunk):
        df_chunk["doubled"] = df_chunk["value"] * 2
        return df_chunk

    df = pd.DataFrame({"value": [1, 2, 3, 4, 5]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, double_column, n_processes=1
    )

    assert len(result) == 5
    assert "doubled" in result.columns
    assert list(result["doubled"]) == [2, 4, 6, 8, 10]


def test_run_pandas_apply_parallel_single_row():
    """Test single-process path with single-row dataframe (edge case)."""
    def add_computed(df_chunk):
        df_chunk["computed"] = df_chunk["x"] + df_chunk["y"]
        return df_chunk

    df = pd.DataFrame({"x": [5], "y": [3]})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, add_computed, n_processes=1
    )

    assert len(result) == 1
    assert result["computed"].iloc[0] == 8


# =============================================================================
# Tests for run_pandas_apply_parallel — multi process
# =============================================================================


def test_run_pandas_apply_parallel_multi_process():
    """Test multi-process path preserves per-row data integrity.

    The implementation shuffles rows before splitting across workers, then
    concatenates results. This verifies every row's relationship survives.
    """
    df = pd.DataFrame({"value": list(range(10))})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, _mp_double_column, n_processes=2
    )

    assert len(result) == 10
    assert "doubled" in result.columns
    assert (result["doubled"] == result["value"] * 2).all()


def test_run_pandas_apply_parallel_multi_process_preserves_columns():
    """Test multi-process path preserves existing columns and adds new ones."""
    df = pd.DataFrame({"a": list(range(6)), "b": list(range(6, 12))})
    result = CRISPRessoMultiProcessing.run_pandas_apply_parallel(
        df, _mp_add_column, n_processes=2
    )

    assert len(result) == 6
    assert (result["new_col"] == "added").all()
    # Original values survived shuffle+split+reassembly
    assert sorted(result["a"]) == list(range(6))
    assert sorted(result["b"]) == list(range(6, 12))


# =============================================================================
# Tests for run_parallel_commands — single process
# =============================================================================


def test_run_parallel_commands_single_process_success():
    """Test single-process path runs all commands successfully."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        ["true", "true", "true"], n_processes=1, descriptor="test"
    )


def test_run_parallel_commands_single_process_failure():
    """Test single-process path raises on command failure."""
    with pytest.raises(Exception, match="was failed"):
        CRISPRessoMultiProcessing.run_parallel_commands(
            ["true", "false", "true"], n_processes=1, descriptor="test"
        )


def test_run_parallel_commands_single_process_continue_on_fail():
    """Test single-process path continues past failures when flag is set."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        ["true", "false", "true"],
        n_processes=1,
        descriptor="test",
        continue_on_fail=True,
    )


def test_run_parallel_commands_single_process_empty_list():
    """Test single-process path with empty command list."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        [], n_processes=1, descriptor="test"
    )


# =============================================================================
# Tests for run_parallel_commands — multi process
# =============================================================================


def test_run_parallel_commands_multi_process_success():
    """Test multi-process path runs all commands successfully."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        ["true", "true", "true", "true"], n_processes=2, descriptor="test"
    )


def test_run_parallel_commands_multi_process_failure():
    """Test multi-process path raises on command failure."""
    with pytest.raises(Exception, match="failed"):
        CRISPRessoMultiProcessing.run_parallel_commands(
            ["true", "false", "true", "true"], n_processes=2, descriptor="test"
        )


def test_run_parallel_commands_multi_process_continue_on_fail():
    """Test multi-process path continues past failures when flag is set."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        ["true", "false", "true", "true"],
        n_processes=2,
        descriptor="test",
        continue_on_fail=True,
    )


def test_run_parallel_commands_multi_process_empty_list():
    """Test multi-process path with empty command list creates pool but succeeds."""
    CRISPRessoMultiProcessing.run_parallel_commands(
        [], n_processes=2, descriptor="test"
    )


# =============================================================================
# Tests for run_crispresso_cmds
# =============================================================================


def test_run_crispresso_cmds_empty_list():
    """Test run_crispresso_cmds with empty list returns immediately."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=[],
        n_processes="1",
        descriptor="test",
    )


def test_run_crispresso_cmds_single_process_success():
    """Test run_crispresso_cmds single-process path."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=["true", "true"],
        n_processes="1",
        descriptor="test",
    )


def test_run_crispresso_cmds_single_process_failure():
    """Test run_crispresso_cmds single-process path raises on failure."""
    with pytest.raises(Exception, match="failed"):
        CRISPRessoMultiProcessing.run_crispresso_cmds(
            crispresso_cmds=["true", "false", "true"],
            n_processes="1",
            descriptor="test",
        )


def test_run_crispresso_cmds_single_process_continue_on_fail():
    """Test run_crispresso_cmds single-process path continues past failures."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=["true", "false", "true"],
        n_processes="1",
        descriptor="test",
        continue_on_fail=True,
    )


def test_run_crispresso_cmds_multi_process_success():
    """Test run_crispresso_cmds multi-process path (imap_unordered)."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=["true", "true", "true"],
        n_processes="2",
        descriptor="test",
    )


def test_run_crispresso_cmds_multi_process_failure():
    """Test run_crispresso_cmds multi-process path raises on failure."""
    with pytest.raises(Exception, match="failed"):
        CRISPRessoMultiProcessing.run_crispresso_cmds(
            crispresso_cmds=["true", "false", "true"],
            n_processes="2",
            descriptor="test",
        )


def test_run_crispresso_cmds_multi_process_continue_on_fail():
    """Test run_crispresso_cmds multi-process path continues past failures."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=["true", "false", "true"],
        n_processes="2",
        descriptor="test",
        continue_on_fail=True,
    )


def test_run_crispresso_cmds_max_processes():
    """Test run_crispresso_cmds with n_processes='max' branch."""
    CRISPRessoMultiProcessing.run_crispresso_cmds(
        crispresso_cmds=["true"],
        n_processes="max",
        descriptor="test",
    )
