"""Focused tests for the shared report layout's web-only navbar additions."""

from pathlib import Path

from jinja2 import Environment, FileSystemLoader


TEMPLATES = Path(__file__).parents[2] / 'CRISPResso2' / 'CRISPRessoReports' / 'templates'


def test_layout_compiles_and_gates_plot_assistant_link():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    template = env.get_template('layout.html')
    source = template.environment.loader.get_source(env, 'layout.html')[0]

    assert 'current_user.is_authenticated and plot_assistant_enabled and plot_assistant_folder_id' in source
    assert 'href="/plot-assistant/{{ plot_assistant_folder_id }}"' in source
    assert 'AI Plot Assistant' in source


def _render_layout(**ctx):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    base = dict(current_user=None, is_default_user=False, is_web=False, C2PRO_INSTALLED=False,
                get_flashed_messages=lambda **kwargs: [])
    base.update(ctx)
    return env.get_template('layout.html').render(**base)


def test_citation_footer_copyright_is_edilytics_on_the_web():
    html = _render_layout(is_web=True, is_default_user=True)

    assert 'If you like CRISPResso2 please support us by citing it' in html
    assert 'Copyright <a href="http://edilytics.com">Edilytics</a>' in html
    assert 'Kendell Clement and Luca Pinello</a>' not in html


def test_citation_footer_copyright_keeps_authors_in_cli_reports():
    html = _render_layout(is_web=False)

    assert 'If you like CRISPResso2 please support us by citing it' in html
    assert 'Copyright <a href="http://pinellolab.org">Kendell Clement and Luca Pinello</a>' in html
    assert 'edilytics.com">Edilytics' not in html
