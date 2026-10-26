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
