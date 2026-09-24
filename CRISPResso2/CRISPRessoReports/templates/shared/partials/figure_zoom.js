/* CRISPRessoReports — hover zoom for wide report figures (figure 2a).
 *
 * Rendered INLINE by shared/partials/figure_zoom.html (reports must stay
 * self-contained: downloaded report zips are opened from file://, so we
 * cannot depend on any external script).
 *
 * The handler is deliberately STATELESS: every pointermove re-measures the
 * figure and the strip and recomputes the lens/pan from scratch. The previous
 * implementation cached geometry at first hover ("hasWidthSet"), which locked
 * in zeros/NaN when the first hover happened inside a hidden tab or before
 * the image loaded — leaving the zoom dead or out of sync for the rest of
 * the page view. Per-event recomputation makes those failure modes
 * impossible: a hidden or unloaded figure simply bails out (computeZoom
 * returns null) and works on the next event once laid out.
 *
 * computeZoom() is unit-tested (numeric cases, via a small node harness)
 * by tests/unit_tests/test_report_figure_zoom.py in this repo.
 */
(function () {
  'use strict';

  /* Pure geometry: given the displayed figure size, the strip size, and the
   * cursor position (figure-relative px), return the lens size/offset and the
   * strip-image translation. Returns null when the geometry is unusable
   * (zero sizes), and {fits: true, ...} when the figure already fits in the
   * strip at strip height (nothing to zoom; centered). */
  function computeZoom(imgW, imgH, stripW, stripH, cursorX) {
    if (!(imgW > 0) || !(imgH > 0) || !(stripW > 0) || !(stripH > 0)) {
      return null;
    }
    var k = stripH / imgH;        /* strip px per figure px */
    var scaledW = imgW * k;       /* figure width when scaled to strip height */
    if (scaledW <= stripW) {
      /* Narrow figure (e.g. a short amplicon quilt): no zoom possible —
       * center it in the strip and let the lens span the whole figure. */
      return {
        fits: true,
        lensW: imgW,
        lensX: 0,
        translateX: Math.round((stripW - scaledW) / 2)
      };
    }
    var lensW = stripW / k;       /* figure px visible through the strip */
    var lensX = cursorX - lensW / 2;
    if (lensX < 0) { lensX = 0; }
    if (lensX > imgW - lensW) { lensX = imgW - lensW; }
    return {
      fits: false,
      lensW: lensW,
      lensX: lensX,
      /* 0 - x avoids producing -0 at the left edge (Object.is(-0, 0) is false) */
      translateX: 0 - Math.round(k * lensX)
    };
  }

  /* Build a quoted attribute-value selector for a uid. Amplicon names come
   * from user FASTA headers and are not guaranteed quote-free; an unescaped
   * quote or backslash would make querySelector throw (killing wiring for
   * every figure not yet wired). Escaping " and \\ suffices inside a
   * quoted CSS string. */
  function selFor(attr, uid) {
    return '[' + attr + '="' + uid.replace(/["\\]/g, '\\$&') + '"]';
  }

  function wireOne(figure) {
    var uid = figure.getAttribute('data-cr-zoom-figure');
    var img = figure.querySelector(selFor('data-cr-zoom-img', uid));
    var lens = figure.querySelector(selFor('data-cr-zoom-lens', uid));
    var strip = figure.ownerDocument.querySelector(selFor('data-cr-zoom-strip', uid));
    var stripImg = strip && strip.querySelector(selFor('data-cr-zoom-strip-img', uid));
    if (!img || !lens || !strip || !stripImg) { return false; }

    function update(e) {
      var z = computeZoom(
        img.offsetWidth, img.offsetHeight,
        strip.offsetWidth, strip.offsetHeight,
        e.clientX - img.getBoundingClientRect().left
      );
      if (!z) { return; }  /* hidden tab / image not laid out yet: retry next event */
      lens.style.width = z.lensW + 'px';
      lens.style.left = z.lensX + 'px';
      stripImg.style.transform = 'translateX(' + z.translateX + 'px)';
    }

    figure.addEventListener('pointermove', update, { passive: true });
    return true;
  }

  function wireAll(doc) {
    var figures = doc.querySelectorAll('[data-cr-zoom-figure]:not([data-cr-zoom-wired])');
    for (var i = 0; i < figures.length; i++) {
      /* stamp only on success: a figure whose elements are missing stays
       * unwired and gets retried by the next include instead of being
       * silently skipped forever. */
      if (wireOne(figures[i])) {
        figures[i].setAttribute('data-cr-zoom-wired', '1');
      }
    }
  }

  /* Browser: wire every zoom figure rendered so far. figure_zoom.html
   * includes this file once per figure, right after its markup, so each
   * inclusion wires any figure not yet wired (idempotent). */
  if (typeof document !== 'undefined' && document.querySelectorAll) {
    wireAll(document);
  }

  /* Test loading: expose the pure function for the node-based geometry
   * tests (TestComputeZoomGeometry in this repo's unit tests). */
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { computeZoom: computeZoom };
  }
})();
