"""Tests for the internal upsetplot module ported from upsetplot v0.9.0."""
import os
import tempfile

import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from CRISPResso2.plots import upsetplot


def _make_sample_data():
    """Create a small boolean-MultiIndex Series like CRISPResso2 produces."""
    index = pd.MultiIndex.from_tuples(
        [
            (True, False, False),
            (False, True, False),
            (True, True, False),
            (False, False, True),
            (True, False, True),
            (False, False, False),
        ],
        names=["3:C->T", "5:A->G", "has_indel"],
    )
    return pd.Series([40, 20, 15, 10, 8, 7], index=index, name="cat_counts")


# ---------------------------------------------------------------------------
# Plot tests
# ---------------------------------------------------------------------------

def test_plot_returns_axes_dict():
    data = _make_sample_data()
    fig = plt.figure(figsize=(10, 6))
    result = upsetplot.plot(
        data,
        fig=fig,
        element_size=None,
        show_counts=True,
        show_percentages="{:.2f}",
        sort_categories_by="-input",
    )
    assert isinstance(result, dict)
    assert "matrix" in result
    assert "intersections" in result
    assert "totals" in result
    assert "shading" in result
    plt.close(fig)


def test_plot_saves_to_pdf():
    data = _make_sample_data()
    fig = plt.figure(figsize=(15, 10))
    upsetplot.plot(
        data,
        fig=fig,
        element_size=None,
        show_counts=True,
        show_percentages="{:.2f}",
        sort_categories_by="-input",
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_upset.pdf")
        fig.savefig(path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0
    plt.close(fig)


def test_plot_saves_to_png():
    data = _make_sample_data()
    fig = plt.figure(figsize=(15, 10))
    upsetplot.plot(
        data,
        fig=fig,
        element_size=None,
        show_counts=True,
        show_percentages="{:.2f}",
        sort_categories_by="-input",
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "test_upset.png")
        fig.savefig(path)
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0
    plt.close(fig)


def test_plot_two_categories():
    """A two-category MultiIndex should produce a valid plot."""
    index = pd.MultiIndex.from_tuples(
        [(True, False), (False, True), (False, False)],
        names=["edit_a", "edit_b"],
    )
    data = pd.Series([30, 50, 20], index=index, name="counts")
    fig = plt.figure(figsize=(10, 6))
    result = upsetplot.plot(data, fig=fig, element_size=None)
    assert "matrix" in result
    plt.close(fig)


def test_plot_all_false_row():
    """Data with a row where all categories are False."""
    index = pd.MultiIndex.from_tuples(
        [
            (True, True),
            (True, False),
            (False, True),
            (False, False),
        ],
        names=["cat_a", "cat_b"],
    )
    data = pd.Series([10, 20, 15, 55], index=index, name="counts")
    fig = plt.figure(figsize=(10, 6))
    result = upsetplot.plot(data, fig=fig, element_size=None)
    assert "matrix" in result
    plt.close(fig)


def test_sort_categories_by_input():
    """sort_categories_by='-input' preserves reverse input order."""
    data = _make_sample_data()
    upset = upsetplot.UpSet(data, sort_categories_by="-input")
    assert list(upset.totals.index) == ["has_indel", "5:A->G", "3:C->T"]


def test_show_counts_and_percentages():
    """Counts and percentages should not raise errors."""
    data = _make_sample_data()
    fig = plt.figure(figsize=(10, 6))
    upsetplot.plot(
        data,
        fig=fig,
        element_size=None,
        show_counts=True,
        show_percentages="{:.1%}",
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# Query tests
# ---------------------------------------------------------------------------

def test_query_basic():
    data = _make_sample_data()
    result = upsetplot.query(data)
    assert result.total == 100
    assert len(result.subset_sizes) == 6
    assert isinstance(result.category_totals, pd.Series)


def test_query_sort_by_cardinality():
    data = _make_sample_data()
    result = upsetplot.query(data, sort_by="cardinality")
    sizes = result.subset_sizes.values
    assert all(sizes[i] >= sizes[i + 1] for i in range(len(sizes) - 1))


def test_query_min_subset_size():
    data = _make_sample_data()
    result = upsetplot.query(data, min_subset_size=10)
    assert all(result.subset_sizes >= 10)


def test_query_non_boolean_index_raises():
    index = pd.MultiIndex.from_tuples(
        [(1, 2), (3, 4)],
        names=["a", "b"],
    )
    data = pd.Series([10, 20], index=index)
    with pytest.raises(ValueError, match="not boolean"):
        upsetplot.query(data)


# ---------------------------------------------------------------------------
# Format conversion tests
# ---------------------------------------------------------------------------

def test_format_conversion_simple_int():
    assert upsetplot._to_new_pos_format("%d") == "{:d}"


def test_format_conversion_float():
    assert upsetplot._to_new_pos_format("%.2f") == "{:.2f}"


def test_format_conversion_passthrough_new_style():
    result = upsetplot._to_new_pos_format("{:.2f}")
    assert ":.2f" in result
