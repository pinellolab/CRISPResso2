"""Render-contract tests for the shared hover-zoom figure macro.

``shared/partials/figure_zoom.html`` is the single hover-zoom implementation
for report figures (legacy plot_2a): it is rendered by the no-Pro
CRISPRessoReports report template AND by CRISPRessoPro's nucleotide-quilt
partial, so its markup and inline behavior are pinned here. The bugs that
motivated the shared implementation were all markup/JS regressions:

- the zoom strip collapsing to 0px wide inside flex-column wrappers
  (missing explicit width)
- one-shot JS state ("hasWidthSet") locking in zero/NaN geometry when the
  first hover happened in a hidden tab or before image load
- Pro and no-Pro renderers drifting apart (different zoom markup, zoom on
  figures where the zoomed copy renders smaller than the plain figure)

The macro must render identically under every jinja environment the report
stack uses: the CLI builder and CRISPRessoPro's hooks (autoescape off) and
the CRISPRessoWEB Flask app (autoescape by file extension). The inline
script is included from figure_zoom.js and must stay unescaped everywhere —
a report zip opened from file:// cannot load external scripts.
"""

import json
import os
import shutil
import subprocess

import pytest
from jinja2 import Environment, FileSystemLoader, select_autoescape
from pathlib import Path

from CRISPResso2 import CRISPRessoReports

TEMPLATES_DIR = os.path.join(os.path.dirname(CRISPRessoReports.__file__), "templates")
JS_PATH = os.path.join(TEMPLATES_DIR, "shared", "partials", "figure_zoom.js")


def _env(flask_like):
    loader = FileSystemLoader(TEMPLATES_DIR)
    if flask_like:
        return Environment(loader=loader, autoescape=select_autoescape(enabled_extensions=("html", "htm", "xml")))
    return Environment(loader=loader)


def _render(flask_like=False, **macro_kwargs):
    kwargs = {"img_src": "out/2a.REF1.Nucleotide_percentage_quilt.png",
              "pdf_href": "out/2a.REF1.Nucleotide_percentage_quilt.pdf",
              "uid": "nucs_REF1"}
    kwargs.update(macro_kwargs)
    env = _env(flask_like)
    template = env.from_string(
        "{% from 'shared/partials/figure_zoom.html' import figure_zoom %}"
        "{{ figure_zoom(img_src, pdf_href, uid, width) }}"
    )
    return template.render(width="95%", **kwargs)


@pytest.mark.parametrize("flask_like", [False, True], ids=["cli/pro-env", "flask-env"])
class TestFigureZoomMarkup:
    def test_strip_has_explicit_width_and_overflow_window(self, flask_like):
        """The strip must carry width:100% + overflow:hidden.

        It is a flex item inside flex-column figure wrappers; without an
        explicit width it shrink-to-fits to 0px (the blank-strip bug) and
        without overflow:hidden the translated image spills instead of
        cropping.
        """
        html = _render(flask_like)
        assert 'data-cr-zoom-strip="nucs_REF1"' in html
        strip_style = html.split('data-cr-zoom-strip="nucs_REF1"', 1)[1].split(">", 1)[0]
        assert "width:100%" in strip_style.replace(" ", "")
        assert "overflow:hidden" in strip_style.replace(" ", "")

    def test_all_five_elements_carry_matching_data_attributes(self, flask_like):
        """Strip / strip-img / figure / lens / img are bound by data-cr-zoom-*."""
        html = _render(flask_like)
        for attr in ("data-cr-zoom-strip=", "data-cr-zoom-strip-img=",
                     "data-cr-zoom-figure=", "data-cr-zoom-lens=", "data-cr-zoom-img="):
            assert attr in html, f"missing {attr} binding"

    def test_pdf_link_wraps_main_image(self, flask_like):
        html = _render(flask_like)
        assert '<a href="out/2a.REF1.Nucleotide_percentage_quilt.pdf">' in html
        assert 'src="out/2a.REF1.Nucleotide_percentage_quilt.png"' in html

    def test_inline_script_is_the_shared_stateless_js(self, flask_like):
        """The figure embeds figure_zoom.js verbatim (reports are self-contained)."""
        html = _render(flask_like)
        js = Path(JS_PATH).read_text()
        assert js in html, "figure_zoom.js must be inlined verbatim by the macro"

    def test_script_is_stateless_no_cached_geometry(self, flask_like):
        """The old one-shot init ('hasWidthSet') bug must not come back."""
        html = _render(flask_like)
        assert "lens.hasWidthSet" not in html
        assert "function computeZoom" in html
        assert "pointermove" in html

    def test_strip_images_defy_bootstrap_img_clamping(self, flask_like):
        """Bootstrap's img{max-width:100%} would clamp and distort the
        full-height strip copy; the macro must opt out explicitly.
        """
        html = _render(flask_like)
        assert html.count("max-width:none") >= 2  # zoom copy + tablet scroll copy

    def test_tablet_scroll_fallback_preserved(self, flask_like):
        """Below lg there is no hover: the strip doubles as a scrollable view."""
        html = _render(flask_like)
        assert 'class="d-lg-none"' in html
        assert "overflow-x:auto" in html.replace(" ", "")

    def test_inline_js_escapes_uid_in_selectors(self, flask_like):
        """Amplicon names come from user FASTA headers — not quote-free.

        The JS reads the uid from the DOM and splices it into quoted
        attribute-value selectors; an unescaped quote or backslash makes
        querySelector throw, which would kill wiring for every figure not
        yet wired on the page. Selectors must route through an escaper.
        """
        html = _render(flask_like)
        assert "function selFor" in html
        # every selector routes through the escaper; the raw-splice form
        # querySelector('[data-cr-zoom-...="' + uid ...) must be gone
        assert html.count("querySelector(selFor(") >= 4
        assert "querySelector('[data-cr-zoom-" not in html


@pytest.mark.parametrize("flask_like", [False, True], ids=["cli/pro-env", "flask-env"])
def test_uid_is_escaped_in_attribute(flask_like):
    """A hostile uid must not break out of the attribute in any environment.

    The CLI/Pro jinja environments run with autoescape OFF, so the macro
    applies |e explicitly (attribute-breaking was possible there before:
    a quoted amplicon name truncated the attribute and broke the JS wiring).
    Under autoescape (the Flask env) |e is idempotent — Markup escapes
    are stable under re-escaping.
    """
    html = _render(flask_like, uid='a"><script>')
    escaped = 'a&#34;&gt;&lt;script&gt;'
    assert escaped in html, 'uid must be HTML-escaped in both environments'
    assert html.count('data-cr-zoom-strip="' + escaped) == 1
    assert html.count('data-cr-zoom-figure="' + escaped) == 1


pytestmark_node = pytest.mark.skipif(
    shutil.which("node") is None, reason="node is not installed"
)


class TestComputeZoomGeometry:
    """Numeric contract for computeZoom — the fix itself — via a node harness.

    computeZoom is pure geometry; these cases pin the zoom math (lens width,
    both clamps, centering) and the two bail-outs (zero geometry -> null,
    fits-in-strip -> centered, no zoom). Expected values are hand-computed
    from the documented semantics:
      k = stripH / imgH;  lensW = stripW / k;  lensX = clamp(cursorX -
      lensW/2, 0, imgW - lensW);  translateX = -round(k * lensX)
    """

    HARNESS = (
        "const { computeZoom } = require(process.argv[1]);"
        "console.log(JSON.stringify(computeZoom.apply(null, JSON.parse(process.argv[2]))));"
    )

    def _compute(self, img_w, img_h, strip_w, strip_h, cursor_x):
        proc = subprocess.run(
            ["node", "-e", self.HARNESS, JS_PATH,
             json.dumps([img_w, img_h, strip_w, strip_h, cursor_x])],
            capture_output=True, text=True, check=True,
        )
        return json.loads(proc.stdout)

    @pytestmark_node
    def test_mid_figure_zoom_unit_scale(self):
        z = self._compute(1000, 100, 200, 100, 500)
        assert z == {"fits": False, "lensW": 200, "lensX": 400, "translateX": -400}

    @pytestmark_node
    def test_mid_figure_zoom_half_scale(self):
        # k = 50/100 = 0.5: lensW = 200/0.5 = 400, lensX = 500-200 = 300,
        # translateX = -round(0.5 * 300) = -150
        z = self._compute(1000, 100, 200, 50, 500)
        assert z == {"fits": False, "lensW": 400, "lensX": 300, "translateX": -150}

    @pytestmark_node
    def test_left_edge_clamps_to_zero(self):
        z = self._compute(1000, 100, 200, 100, 0)
        assert z == {"fits": False, "lensW": 200, "lensX": 0, "translateX": 0}

    @pytestmark_node
    def test_right_edge_clamps_to_imgW_minus_lensW(self):
        z = self._compute(1000, 100, 200, 100, 1000)
        assert z == {"fits": False, "lensW": 200, "lensX": 800, "translateX": -800}

    @pytestmark_node
    def test_zero_geometry_returns_null(self):
        # hidden tab / not laid out yet: the caller retries on the next event
        assert self._compute(1000, 100, 0, 100, 500) is None
        assert self._compute(0, 100, 200, 100, 500) is None

    @pytestmark_node
    def test_narrow_figure_fits_and_centers(self):
        # scaledW = 200 * (100/100) = 200 <= stripW 400: nothing to zoom,
        # the scaled figure is centered instead
        z = self._compute(200, 100, 400, 100, 999)
        assert z == {"fits": True, "lensW": 200, "lensX": 0, "translateX": 100}


class TestReportZoomWiring:
    """report.html must wire plot_2a locs through the shared macro (and must
    not resurrect the removed per-amplicon zoom ids/JS it used to inline).
    """

    def test_plot_2a_locs_render_shared_zoom(self):
        from CRISPResso2.CRISPRessoReports.jinja_partials import (
            generate_render_partial, render_partial,
        )

        env = _env(False)

        def custom_partial_render(name, **data):
            template = env.get_template(name)
            data.update(
                render_partial=generate_render_partial(custom_partial_render),
                is_default_user=False, is_web=False, C2PRO_INSTALLED=False,
            )
            return template.render(**data)

        report_data = {
            "amplicons": ["REF1"],
            "figures": {
                "htmls": {"REF1": {}},
                "locs": {"REF1": {"plot_2a": "2a.REF1.Nucleotide_percentage_quilt"}},
                "titles": {"REF1": {"plot_2a": "Figure 2a"}},
                "captions": {"REF1": {"plot_2a": "cap"}},
                "datas": {"REF1": {"plot_2a": []}},
                "names": {"REF1": ["plot_2a"]},
                "sgRNA_based_names": {"REF1": {}},
            },
            "run_data": {"running_info": {"args": {}},
                         "results": {"refs": {"REF1": {}}, "ref_names": ["REF1"],
                                     "general_plots": {}}},
            "report_display_name": "test", "crispresso_data_path": "",
            "nuc_quilt_names": [], "allele_table_names": [],
        }
        html = render_partial("report.html", custom_partial_render,
                              report_data=report_data, C2PRO_INSTALLED=False)
        assert 'data-cr-zoom-figure="nucs_REF1"' in html
        assert 'src="2a.REF1.Nucleotide_percentage_quilt.png"' in html
        assert 'href="2a.REF1.Nucleotide_percentage_quilt.pdf"' in html
        assert html.count("function computeZoom") == 1
        # the retired per-amplicon zoom machinery must stay gone
        for gone in ("tozoom_nucs_REF1", "zoomview_nucs_REF1",
                     "zoomlens_nucs_REF1", "updateZoom"):
            assert gone not in html
