"""Minimal UpSet plot implementation for CRISPResso2.

Ported from UpSetPlot v0.9.0 (https://github.com/jnothman/UpSetPlot)
by Joel Nothman and contributors.

This code was ported into CRISPResso2 because the upstream upsetplot package
is no longer maintained and is incompatible with newer versions of pandas (>2),
matplotlib (>=3.8), and NumPy (>=1.24).

Only the subset of functionality used by CRISPResso2 is included here:
- The ``plot()`` convenience function
- The ``UpSet`` class (without catplot, stacked bars, style_subsets,
  style_categories)
- Supporting data reformatting and query logic

Original license (BSD):

    Copyright (c) 2018-2024 Joel Nothman.
    All rights reserved.


    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    a. Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
    b. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
    c. The names of the contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.


    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.
"""

import re
import typing
import warnings

import matplotlib
import matplotlib.gridspec
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib import pyplot as plt

# matplotlib.tight_layout.get_renderer was removed in matplotlib >=3.6.
try:
    from matplotlib.tight_layout import get_renderer
    _RENDERER_IMPORTED = True
except ImportError:
    _RENDERER_IMPORTED = False


# ---------------------------------------------------------------------------
# Format-string conversion utilities  (from upsetplot/util.py)
# ---------------------------------------------------------------------------

_ODD_REPEAT_PATTERN = r"((?<!{c}){c}({c}{c})*(?!{c}))"
_SPECIFIER_PATTERN = r"[^(]*[diouxXeEfFgGcs]"


def _to_new_pos_format(fmt: str) -> str:
    """Convert old-style positional ``%``-formatting to new-style ``{}``."""
    odd_perc_pattern = _ODD_REPEAT_PATTERN.format(c="%")
    fmt = fmt.replace("{", "{{").replace("}", "}}")

    even_perc_pattern = r"(?<!%)(%%)+" r"(?!%)s"
    if re.search(even_perc_pattern, fmt):
        raise TypeError(
            "not all arguments converted during string formatting"
        )
    fmt = re.sub(
        rf"%({_SPECIFIER_PATTERN})",
        lambda m: f"{{:{m.group(1)}}}",
        fmt,
    )
    if re.search(odd_perc_pattern, fmt):
        raise ValueError("incomplete format")
    return fmt.replace("%%", "%")


# ---------------------------------------------------------------------------
# Data aggregation / query helpers  (from upsetplot/reformat.py)
# ---------------------------------------------------------------------------

class QueryResult:
    """Container for reformatted data and aggregates."""

    def __init__(self, data, subset_sizes, category_totals, total):
        self.data = data
        self.subset_sizes = subset_sizes
        self.category_totals = category_totals
        self.total = total


def _aggregate_data(df, subset_size, sum_over):
    _SUBSET_SIZE_VALUES = ["auto", "count", "sum"]
    if subset_size not in _SUBSET_SIZE_VALUES:
        raise ValueError(
            f"subset_size should be one of {_SUBSET_SIZE_VALUES}. "
            f"Got {subset_size!r}"
        )
    if df.ndim == 1:
        input_name = df.name
        df = pd.DataFrame({"_value": df})
        if subset_size == "auto" and not df.index.is_unique:
            raise ValueError(
                'subset_size="auto" cannot be used for a '
                "Series with non-unique groups."
            )
        if sum_over is not None:
            raise ValueError(
                "sum_over is not applicable when the input is a Series"
            )
        sum_over = False if subset_size == "count" else "_value"
    elif sum_over is False:
        raise ValueError("Unsupported value for sum_over: False")
    elif subset_size == "auto" and sum_over is None:
        sum_over = False
    elif subset_size == "count":
        if sum_over is not None:
            raise ValueError(
                "sum_over cannot be set if subset_size=%r" % subset_size
            )
        sum_over = False
    elif subset_size == "sum" and sum_over is None:
        raise ValueError(
            "sum_over should be a field name if "
            'subset_size="sum" and a DataFrame is provided.'
        )

    gb = df.groupby(level=list(range(df.index.nlevels)), sort=False)
    if sum_over is False:
        aggregated = gb.size()
        aggregated.name = "size"
    elif hasattr(sum_over, "lower"):
        aggregated = gb[sum_over].sum()
    else:
        raise ValueError("Unsupported value for sum_over: %r" % sum_over)

    if aggregated.name == "_value":
        aggregated.name = input_name

    return df, aggregated


def _check_index(df):
    if not all({True, False} >= set(level) for level in df.index.levels):
        raise ValueError(
            "The DataFrame has values in its index that are not boolean"
        )
    df = df.copy(deep=False)
    kw = {
        "levels": [x.astype(bool) for x in df.index.levels],
        "codes": df.index.codes,
        "names": df.index.names,
    }
    df.index = pd.MultiIndex(**kw)
    return df


def _scalar_to_list(val):
    if not isinstance(val, (typing.Sequence, set)) or isinstance(val, str):
        val = [val]
    return val


def _check_percent(value, agg):
    if not isinstance(value, str):
        return value
    try:
        if value.endswith("%") and 0 <= float(value[:-1]) <= 100:
            return float(value[:-1]) / 100 * agg.sum()
    except ValueError:
        pass
    raise ValueError(
        "String value must be formatted as percentage between 0 and 100. "
        f"Got {value}"
    )


def _get_subset_mask(
    agg,
    min_subset_size,
    max_subset_size,
    max_subset_rank,
    min_degree,
    max_degree,
    present,
    absent,
):
    min_subset_size = _check_percent(min_subset_size, agg)
    max_subset_size = _check_percent(max_subset_size, agg)
    subset_mask = True
    if min_subset_size is not None:
        subset_mask = np.logical_and(subset_mask, agg >= min_subset_size)
    if max_subset_size is not None:
        subset_mask = np.logical_and(subset_mask, agg <= max_subset_size)
    if max_subset_rank is not None:
        subset_mask = np.logical_and(
            subset_mask,
            agg.rank(method="min", ascending=False) <= max_subset_rank,
        )
    if (min_degree is not None and min_degree >= 0) or max_degree is not None:
        degree = agg.index.to_frame().sum(axis=1)
        if min_degree is not None:
            subset_mask = np.logical_and(subset_mask, degree >= min_degree)
        if max_degree is not None:
            subset_mask = np.logical_and(subset_mask, degree <= max_degree)
    if present is not None:
        for col in _scalar_to_list(present):
            subset_mask = np.logical_and(
                subset_mask, agg.index.get_level_values(col).values
            )
    if absent is not None:
        for col in _scalar_to_list(absent):
            exclude_mask = np.logical_not(
                agg.index.get_level_values(col).values
            )
            subset_mask = np.logical_and(subset_mask, exclude_mask)
    return subset_mask


def _filter_subsets(
    df,
    agg,
    min_subset_size,
    max_subset_size,
    max_subset_rank,
    min_degree,
    max_degree,
    present,
    absent,
):
    subset_mask = _get_subset_mask(
        agg,
        min_subset_size=min_subset_size,
        max_subset_size=max_subset_size,
        max_subset_rank=max_subset_rank,
        min_degree=min_degree,
        max_degree=max_degree,
        present=present,
        absent=absent,
    )
    if subset_mask is True:
        return df, agg
    agg = agg[subset_mask]
    df = df[df.index.isin(agg.index)]
    return df, agg


def query(
    data,
    sort_by="degree",
    sort_categories_by="cardinality",
    subset_size="auto",
    sum_over=None,
    min_subset_size=None,
    max_subset_size=None,
    max_subset_rank=None,
    min_degree=None,
    max_degree=None,
    present=None,
    absent=None,
    include_empty_subsets=False,
):
    """Transform and filter a categorised dataset."""
    data, agg = _aggregate_data(data, subset_size, sum_over)
    data = _check_index(data)
    grand_total = agg.sum()
    category_totals = [
        agg[agg.index.get_level_values(name).values.astype(bool)].sum()
        for name in agg.index.names
    ]
    category_totals = pd.Series(category_totals, index=agg.index.names)

    if include_empty_subsets:
        nlevels = len(agg.index.levels)
        if nlevels > 10:
            raise ValueError(
                "include_empty_subsets is supported for at most 10 categories"
            )
        new_agg = pd.Series(
            0,
            index=pd.MultiIndex.from_product(
                [[False, True]] * nlevels, names=agg.index.names
            ),
            dtype=agg.dtype,
            name=agg.name,
        )
        new_agg.update(agg)
        agg = new_agg

    data, agg = _filter_subsets(
        data,
        agg,
        min_subset_size=min_subset_size,
        max_subset_size=max_subset_size,
        max_subset_rank=max_subset_rank,
        min_degree=min_degree,
        max_degree=max_degree,
        present=present,
        absent=absent,
    )

    # Sort categories
    if sort_categories_by in ("cardinality", "-cardinality"):
        category_totals.sort_values(
            ascending=sort_categories_by[:1] == "-", inplace=True
        )
    elif sort_categories_by == "-input":
        category_totals = category_totals[::-1]
    elif sort_categories_by in (None, "input"):
        pass
    else:
        raise ValueError(
            "Unknown sort_categories_by: %r" % sort_categories_by
        )
    data = data.reorder_levels(category_totals.index.values)
    agg = agg.reorder_levels(category_totals.index.values)

    # Sort subsets
    if sort_by in ("cardinality", "-cardinality"):
        agg = agg.sort_values(ascending=sort_by[:1] == "-")
    elif sort_by in ("degree", "-degree"):
        index_tuples = sorted(
            agg.index,
            key=lambda x: (sum(x),) + tuple(reversed(x)),
            reverse=sort_by[:1] == "-",
        )
        agg = agg.reindex(
            pd.MultiIndex.from_tuples(index_tuples, names=agg.index.names)
        )
    elif sort_by == "-input":
        agg = agg[::-1]
    elif sort_by in (None, "input"):
        pass
    else:
        raise ValueError("Unknown sort_by: %r" % sort_by)

    return QueryResult(
        data=data,
        subset_sizes=agg,
        category_totals=category_totals,
        total=grand_total,
    )


# ---------------------------------------------------------------------------
# Plotting helpers  (from upsetplot/plotting.py)
# ---------------------------------------------------------------------------

def _process_data(
    df,
    *,
    sort_by,
    sort_categories_by,
    subset_size,
    sum_over,
    min_subset_size=None,
    max_subset_size=None,
    max_subset_rank=None,
    min_degree=None,
    max_degree=None,
    reverse=False,
    include_empty_subsets=False,
):
    results = query(
        df,
        sort_by=sort_by,
        sort_categories_by=sort_categories_by,
        subset_size=subset_size,
        sum_over=sum_over,
        min_subset_size=min_subset_size,
        max_subset_size=max_subset_size,
        max_subset_rank=max_subset_rank,
        min_degree=min_degree,
        max_degree=max_degree,
        include_empty_subsets=include_empty_subsets,
    )

    df = results.data
    agg = results.subset_sizes

    def _pack_binary(X):
        X = pd.DataFrame(X)
        # Use Python object dtype when arbitrary-precision integers are needed.
        dtype = object if X.shape[1] > 62 else np.uint64
        out = pd.Series(0, index=X.index, dtype=dtype)
        for _, col in X.items():
            out *= 2
            out += col
        return out

    df_packed = _pack_binary(df.index.to_frame())
    data_packed = _pack_binary(agg.index.to_frame())
    df["_bin"] = pd.Series(df_packed).map(
        pd.Series(
            np.arange(len(data_packed))[:: -1 if reverse else 1],
            index=data_packed,
        )
    )
    if reverse:
        agg = agg[::-1]

    return results.total, df, agg, results.category_totals


def _multiply_alpha(c, mult):
    r, g, b, a = colors.to_rgba(c)
    a *= mult
    return colors.to_hex((r, g, b, a), keep_alpha=True)


class _Transposed:
    """Wrap an object to transpose horizontal/vertical plotting operations."""

    def __init__(self, obj):
        self.__obj = obj

    def __getattr__(self, key):
        return getattr(self.__obj, self._NAME_TRANSPOSE.get(key, key))

    def __call__(self, *args, **kwargs):
        return self.__obj(
            *args,
            **{self._NAME_TRANSPOSE.get(k, k): v for k, v in kwargs.items()},
        )

    _NAME_TRANSPOSE = {
        "align_xlabels": "align_ylabels",
        "align_ylabels": "align_xlabels",
        "bar": "barh",
        "barh": "bar",
        "bottom": "left",
        "get_figheight": "get_figwidth",
        "get_figwidth": "get_figheight",
        "get_xlim": "get_ylim",
        "get_ylim": "get_xlim",
        "height": "width",
        "hlines": "vlines",
        "hspace": "wspace",
        "left": "bottom",
        "right": "top",
        "set_autoscalex_on": "set_autoscaley_on",
        "set_autoscaley_on": "set_autoscalex_on",
        "set_figheight": "set_figwidth",
        "set_figwidth": "set_figheight",
        "set_xlabel": "set_ylabel",
        "set_xlim": "set_ylim",
        "set_ylabel": "set_xlabel",
        "set_ylim": "set_xlim",
        "sharex": "sharey",
        "sharey": "sharex",
        "top": "right",
        "vlines": "hlines",
        "width": "height",
        "wspace": "hspace",
        "xaxis": "yaxis",
        "yaxis": "xaxis",
    }


def _transpose(obj):
    if isinstance(obj, str):
        return _Transposed._NAME_TRANSPOSE.get(obj, obj)
    return _Transposed(obj)


def _identity(obj):
    return obj


# ---------------------------------------------------------------------------
# UpSet class  (from upsetplot/plotting.py, trimmed to essentials)
# ---------------------------------------------------------------------------

class UpSet:
    """Manage the data and drawing for an UpSet plot.

    This is a minimal version that supports only the features used by
    CRISPResso2: horizontal orientation, bar + matrix + totals, counts and
    percentage labels.

    Parameters
    ----------
    data : pandas.Series or pandas.DataFrame
        Elements associated with categories (a DataFrame), or the size of
        each subset of categories (a Series).  Should have MultiIndex where
        each level is binary, corresponding to category membership.
    orientation : {'horizontal', 'vertical'}
    sort_by : str
    sort_categories_by : str
    subset_size : str
    sum_over : str or None
    min_subset_size, max_subset_size, max_subset_rank : optional
    min_degree, max_degree : optional
    facecolor : str
    other_dots_color : matplotlib color or float
    shading_color : matplotlib color or float
    with_lines : bool
    element_size : float or None
    intersection_plot_elements : int
    totals_plot_elements : int
    show_counts : bool or str
    show_percentages : bool or str
    include_empty_subsets : bool

    """

    _default_figsize = (10, 6)
    DPI = 100

    def __init__(
        self,
        data,
        orientation="horizontal",
        sort_by="degree",
        sort_categories_by="cardinality",
        subset_size="auto",
        sum_over=None,
        min_subset_size=None,
        max_subset_size=None,
        max_subset_rank=None,
        min_degree=None,
        max_degree=None,
        facecolor="auto",
        other_dots_color=0.18,
        shading_color=0.05,
        with_lines=True,
        element_size=32,
        intersection_plot_elements=6,
        totals_plot_elements=2,
        show_counts="",
        show_percentages=False,
        include_empty_subsets=False,
    ):
        self._horizontal = orientation == "horizontal"
        self._reorient = _identity if self._horizontal else _transpose
        if facecolor == "auto":
            bgcolor = matplotlib.rcParams.get("axes.facecolor", "white")
            r, g, b, a = colors.to_rgba(bgcolor)
            lightness = colors.rgb_to_hsv((r, g, b))[-1] * a
            facecolor = "black" if lightness >= 0.5 else "white"
        self._facecolor = facecolor
        self._shading_color = (
            _multiply_alpha(facecolor, shading_color)
            if isinstance(shading_color, float)
            else shading_color
        )
        self._other_dots_color = (
            _multiply_alpha(facecolor, other_dots_color)
            if isinstance(other_dots_color, float)
            else other_dots_color
        )
        self._with_lines = with_lines
        self._element_size = element_size
        self._totals_plot_elements = totals_plot_elements
        self._subset_plots = [
            {
                "type": "default",
                "id": "intersections",
                "elements": intersection_plot_elements,
            }
        ]
        if not intersection_plot_elements:
            self._subset_plots.pop()
        self._show_counts = show_counts
        self._show_percentages = show_percentages

        (self.total, self._df, self.intersections, self.totals) = _process_data(
            data,
            sort_by=sort_by,
            sort_categories_by=sort_categories_by,
            subset_size=subset_size,
            sum_over=sum_over,
            min_subset_size=min_subset_size,
            max_subset_size=max_subset_size,
            max_subset_rank=max_subset_rank,
            min_degree=min_degree,
            max_degree=max_degree,
            reverse=not self._horizontal,
            include_empty_subsets=include_empty_subsets,
        )
        self.subset_styles = [
            {"facecolor": facecolor} for _ in range(len(self.intersections))
        ]

    def _swapaxes(self, x, y):
        if self._horizontal:
            return x, y
        return y, x

    # ------------------------------------------------------------------
    # Bar helpers
    # ------------------------------------------------------------------

    def _plot_bars(self, ax, data, title, colors=None, use_labels=False):
        ax = self._reorient(ax)
        ax.set_autoscalex_on(False)
        data_df = pd.DataFrame(data)
        if self._horizontal:
            data_df = data_df.loc[:, ::-1]

        if callable(colors):
            colors = colors(range(data_df.shape[1]))
        elif isinstance(colors, (str, type(None))):
            colors = [colors] * len(data_df)

        if self._horizontal:
            colors = list(reversed(colors))

        x = np.arange(len(data_df))
        cum_y = None
        all_rects = []
        for (name, y), color in zip(data_df.items(), colors):
            rects = ax.bar(
                x,
                y,
                0.5,
                cum_y,
                color=color,
                zorder=10,
                label=name if use_labels else None,
                align="center",
            )
            cum_y = y if cum_y is None else cum_y + y
            all_rects.extend(rects)

        self._label_sizes(ax, rects, "top" if self._horizontal else "right")

        ax.xaxis.set_visible(False)
        for x in ["top", "bottom", "right"]:
            ax.spines[self._reorient(x)].set_visible(False)

        tick_axis = ax.yaxis
        tick_axis.grid(True)
        ax.set_ylabel(title)
        return all_rects

    def _label_sizes(self, ax, rects, where):
        if not self._show_counts and not self._show_percentages:
            return
        if self._show_counts is True:
            count_fmt = "{:.0f}"
        else:
            count_fmt = self._show_counts
            if count_fmt and "{" not in count_fmt:
                count_fmt = _to_new_pos_format(count_fmt)

        pct_fmt = (
            "{:.1%}"
            if self._show_percentages is True
            else self._show_percentages
        )

        if count_fmt and pct_fmt:
            if where == "top":
                fmt = f"{count_fmt}\n({pct_fmt})"
            else:
                fmt = f"{count_fmt} ({pct_fmt})"

            def make_args(val):
                return val, val / self.total
        elif count_fmt:
            fmt = count_fmt

            def make_args(val):
                return (val,)
        else:
            fmt = pct_fmt

            def make_args(val):
                return (val / self.total,)

        if where == "right":
            margin = 0.01 * abs(np.diff(ax.get_xlim()))
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(
                    float(np.ravel(width + margin)[0]),
                    float(np.ravel(rect.get_y() + rect.get_height() * 0.5)[0]),
                    fmt.format(*make_args(width)),
                    ha="left",
                    va="center",
                )
        elif where == "left":
            margin = 0.01 * abs(np.diff(ax.get_xlim()))
            for rect in rects:
                width = rect.get_width() + rect.get_x()
                ax.text(
                    float(np.ravel(width + margin)[0]),
                    float(np.ravel(rect.get_y() + rect.get_height() * 0.5)[0]),
                    fmt.format(*make_args(width)),
                    ha="right",
                    va="center",
                )
        elif where == "top":
            margin = 0.01 * abs(np.diff(ax.get_ylim()))
            for rect in rects:
                height = rect.get_height() + rect.get_y()
                ax.text(
                    float(np.ravel(rect.get_x() + rect.get_width() * 0.5)[0]),
                    float(np.ravel(height + margin)[0]),
                    fmt.format(*make_args(height)),
                    ha="center",
                    va="bottom",
                )
        else:
            raise NotImplementedError("unhandled where: %r" % where)

    # ------------------------------------------------------------------
    # Grid layout
    # ------------------------------------------------------------------

    def make_grid(self, fig=None):
        """Get a SubplotSpec for each Axes, accounting for label text width."""
        n_cats = len(self.totals)
        n_inters = len(self.intersections)

        if fig is None:
            fig = plt.gcf()

        text_kw = {"size": matplotlib.rcParams["xtick.labelsize"]}
        t = fig.text(
            0,
            0,
            "\n".join(str(label) + "x" for label in self.totals.index.values),
            **text_kw,
        )
        window_extent_args = {}
        if _RENDERER_IMPORTED:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                window_extent_args["renderer"] = get_renderer(fig)
        textw = t.get_window_extent(**window_extent_args).width
        t.remove()

        window_extent_args = {}
        if _RENDERER_IMPORTED:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                window_extent_args["renderer"] = get_renderer(fig)
        figw = self._reorient(
            fig.get_window_extent(**window_extent_args)
        ).width

        sizes = np.asarray([p["elements"] for p in self._subset_plots])
        fig = self._reorient(fig)

        non_text_nelems = len(self.intersections) + self._totals_plot_elements
        if self._element_size is None:
            colw = (figw - textw) / non_text_nelems
        else:
            render_ratio = figw / fig.get_figwidth()
            colw = self._element_size / 72 * render_ratio
            figw = colw * (non_text_nelems + np.ceil(textw / colw) + 1)
            fig.set_figwidth(figw / render_ratio)
            fig.set_figheight((colw * (n_cats + sizes.sum())) / render_ratio)

        text_nelems = int(np.ceil(figw / colw - non_text_nelems))

        GS = self._reorient(matplotlib.gridspec.GridSpec)
        gridspec = GS(
            *self._swapaxes(
                n_cats + (sizes.sum() or 0),
                n_inters + text_nelems + self._totals_plot_elements,
            ),
            hspace=1,
        )
        if self._horizontal:
            out = {
                "matrix": gridspec[-n_cats:, -n_inters:],
                "shading": gridspec[-n_cats:, :],
                "totals": (
                    None
                    if self._totals_plot_elements == 0
                    else gridspec[-n_cats:, : self._totals_plot_elements]
                ),
                "gs": gridspec,
            }
            cumsizes = np.cumsum(sizes[::-1])
            for start, stop, p in zip(
                np.hstack([[0], cumsizes]),
                cumsizes,
                self._subset_plots[::-1],
            ):
                out[p["id"]] = gridspec[start:stop, -n_inters:]
        else:
            out = {
                "matrix": gridspec[-n_inters:, :n_cats],
                "shading": gridspec[:, :n_cats],
                "totals": (
                    None
                    if self._totals_plot_elements == 0
                    else gridspec[: self._totals_plot_elements, :n_cats]
                ),
                "gs": gridspec,
            }
            cumsizes = np.cumsum(sizes)
            for start, stop, p in zip(
                np.hstack([[0], cumsizes]),
                cumsizes,
                self._subset_plots,
            ):
                out[p["id"]] = gridspec[
                    -n_inters:, start + n_cats : stop + n_cats
                ]
        return out

    # ------------------------------------------------------------------
    # Matrix / intersection dots
    # ------------------------------------------------------------------

    def plot_matrix(self, ax):
        """Plot the matrix of intersection indicators onto *ax*."""
        ax = self._reorient(ax)
        data = self.intersections
        n_cats = data.index.nlevels

        inclusion = data.index.to_frame().values

        styles = [
            [
                self.subset_styles[i]
                if inclusion[i, j]
                else {
                    "facecolor": self._other_dots_color,
                    "linewidth": 0,
                }
                for j in range(n_cats)
            ]
            for i in range(len(data))
        ]
        styles = sum(styles, [])  # flatten
        style_columns = {
            "facecolor": "facecolors",
            "edgecolor": "edgecolors",
            "linewidth": "linewidths",
            "linestyle": "linestyles",
            "hatch": "hatch",
        }
        styles = (
            pd.DataFrame(styles)
            .reindex(columns=style_columns.keys())
            .astype(
                {
                    "facecolor": "O",
                    "edgecolor": "O",
                    "linewidth": float,
                    "linestyle": "O",
                    "hatch": "O",
                }
            )
        )
        styles["linewidth"] = styles["linewidth"].fillna(1)
        styles["facecolor"] = styles["facecolor"].fillna(self._facecolor)
        styles["edgecolor"] = styles["edgecolor"].fillna(styles["facecolor"])
        styles["linestyle"] = styles["linestyle"].fillna("solid")
        del styles["hatch"]

        x = np.repeat(np.arange(len(data)), n_cats)
        y = np.tile(np.arange(n_cats), len(data))

        if self._element_size is not None:
            s = (self._element_size * 0.35) ** 2
        else:
            s = 200
        ax.scatter(
            *self._swapaxes(x, y),
            s=s,
            zorder=10,
            **styles.rename(columns=style_columns),
        )

        if self._with_lines:
            idx = np.flatnonzero(inclusion)
            line_data = (
                pd.Series(y[idx], index=x[idx])
                .groupby(level=0)
                .aggregate(["min", "max"])
            )
            line_colors = pd.Series(
                [
                    style.get(
                        "edgecolor",
                        style.get("facecolor", self._facecolor),
                    )
                    for style in self.subset_styles
                ],
                name="color",
            )
            line_data = line_data.join(line_colors)
            ax.vlines(
                line_data.index.values,
                line_data["min"],
                line_data["max"],
                lw=2,
                colors=line_data["color"],
                zorder=5,
            )

        tick_axis = ax.yaxis
        tick_axis.set_ticks(np.arange(n_cats))
        tick_axis.set_ticklabels(
            data.index.names,
            rotation=0 if self._horizontal else -90,
        )
        ax.xaxis.set_visible(False)
        ax.tick_params(axis="both", which="both", length=0)
        if not self._horizontal:
            ax.yaxis.set_ticks_position("top")
        ax.set_frame_on(False)
        ax.set_xlim(-0.5, x[-1] + 0.5, auto=False)
        ax.grid(False)

    # ------------------------------------------------------------------
    # Intersection-size bars
    # ------------------------------------------------------------------

    def plot_intersections(self, ax):
        """Plot bars indicating intersection size."""
        rects = self._plot_bars(
            ax,
            self.intersections,
            title="Intersection size",
            colors=self._facecolor,
        )
        for style, rect in zip(self.subset_styles, rects):
            style = style.copy()
            style.setdefault(
                "edgecolor", style.get("facecolor", self._facecolor)
            )
            for attr, val in style.items():
                getattr(rect, "set_" + attr)(val)

    # ------------------------------------------------------------------
    # Totals bars
    # ------------------------------------------------------------------

    def plot_totals(self, ax):
        """Plot bars indicating total set size."""
        orig_ax = ax
        ax = self._reorient(ax)
        rects = ax.barh(
            np.arange(len(self.totals.index.values)),
            self.totals,
            0.5,
            color=self._facecolor,
            align="center",
        )
        self._label_sizes(
            ax, rects, "left" if self._horizontal else "top"
        )

        max_total = self.totals.max()
        if self._horizontal:
            orig_ax.set_xlim(max_total, 0)
        for x in ["top", "left", "right"]:
            ax.spines[self._reorient(x)].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.xaxis.grid(True)
        ax.yaxis.grid(False)
        ax.patch.set_visible(False)

    # ------------------------------------------------------------------
    # Shading
    # ------------------------------------------------------------------

    def plot_shading(self, ax):
        for i in range(len(self.totals)):
            default_shading = (
                self._shading_color
                if i % 2 == 0
                else (0.0, 0.0, 0.0, 0.0)
            )

            lw = 0
            lw_padding = lw / (self._default_figsize[0] * self.DPI)
            start_x = lw_padding
            end_x = 1 - lw_padding * 3

            rect = plt.Rectangle(
                self._swapaxes(start_x, i - 0.4),
                *self._swapaxes(end_x, 0.8),
                facecolor=default_shading,
                lw=lw,
                zorder=0,
            )
            ax.add_patch(rect)

        ax.set_frame_on(False)
        ax.tick_params(
            axis="both",
            which="both",
            left=False,
            right=False,
            bottom=False,
            top=False,
            labelbottom=False,
            labelleft=False,
        )
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    # ------------------------------------------------------------------
    # Main plot entry point
    # ------------------------------------------------------------------

    def plot(self, fig=None):
        """Draw all parts of the plot onto *fig* or a new figure.

        Returns
        -------
        subplots : dict of matplotlib.axes.Axes
            Keys are 'matrix', 'intersections', 'totals', 'shading'.

        """
        if fig is None:
            fig = plt.figure(figsize=self._default_figsize)
        specs = self.make_grid(fig)
        shading_ax = fig.add_subplot(specs["shading"])
        self.plot_shading(shading_ax)
        matrix_ax = self._reorient(fig.add_subplot)(
            specs["matrix"], sharey=shading_ax
        )
        self.plot_matrix(matrix_ax)
        if specs["totals"] is None:
            totals_ax = None
        else:
            totals_ax = self._reorient(fig.add_subplot)(
                specs["totals"], sharey=matrix_ax
            )
            self.plot_totals(totals_ax)
        out = {
            "matrix": matrix_ax,
            "shading": shading_ax,
            "totals": totals_ax,
        }

        for p in self._subset_plots:
            ax = self._reorient(fig.add_subplot)(
                specs[p["id"]], sharex=matrix_ax
            )
            if p["type"] == "default":
                self.plot_intersections(ax)
            else:
                raise ValueError(
                    "Unknown subset plot type: %r" % p["type"]
                )
            out[p["id"]] = ax

        self._reorient(fig).align_ylabels(
            [out[p["id"]] for p in self._subset_plots]
        )
        return out


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def plot(data, fig=None, **kwargs):
    """Make an UpSet plot of *data* on *fig*.

    Parameters
    ----------
    data : pandas.Series or pandas.DataFrame
        Values for each set to plot.  Should have multi-index where each
        level is binary, corresponding to set membership.
    fig : matplotlib.figure.Figure, optional
        Defaults to a new figure.
    kwargs
        Other arguments for :class:`UpSet`.

    Returns
    -------
    subplots : dict of matplotlib.axes.Axes

    """
    return UpSet(data, **kwargs).plot(fig)
