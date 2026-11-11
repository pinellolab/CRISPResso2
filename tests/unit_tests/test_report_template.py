"""Focused tests for the shared report layout's web-only navbar additions."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from jinja2 import ChoiceLoader, DictLoader, Environment, FileSystemLoader


TEMPLATES = Path(__file__).parents[2] / 'CRISPResso2' / 'CRISPRessoReports' / 'templates'


def test_layout_compiles_and_gates_plot_assistant_link():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    template = env.get_template('layout.html')
    source = template.environment.loader.get_source(env, 'layout.html')[0]

    assert 'current_user.is_authenticated and plot_assistant_enabled' in source
    assert 'plot_assistant_url or' in source
    assert 'AI Plot Assistant' in source


@pytest.mark.parametrize('child', [False, True])
@pytest.mark.parametrize('gate', ['enabled', 'disabled', 'default-user', 'anonymous', 'cli'])
def test_layout_assistant_targets_and_visibility(child, gate):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    context = {
        'is_web': gate != 'cli',
        'is_default_user': gate == 'default-user',
        'current_user': SimpleNamespace(is_authenticated=gate != 'anonymous', email='test@example.com', role='User'),
        'config': {'BANNER_TEXT': 'CRISPResso'},
        'get_flashed_messages': lambda **kwargs: [],
        'url_for': lambda endpoint: '/' + endpoint,
        'plot_assistant_enabled': gate != 'disabled',
        'plot_assistant_folder_id': 'parent',
    }
    expected = '/plot-assistant/parent'
    if child:
        expected += '/children/sample'
        context['plot_assistant_url'] = expected
    # Without an explicit URL, older C2Web versions still get a standalone link.
    html = env.get_template('layout.html').render(**context)
    if gate == 'enabled':
        assert f'href="{expected}"' in html
        assert html.count('AI Plot Assistant') == 1
        if child:
            assert 'href="/plot-assistant/parent"' not in html
    else:
        assert 'AI Plot Assistant' not in html


def test_core_report_header_uses_feature_flag_and_explicit_target():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    source = env.loader.get_source(env, 'report.html')[0]
    assert "report_data.get('plot_assistant_enabled')" in source
    assert 'href="{{ report_data[\'plot_assistant_url\'] }}"' in source
    assert 'href="/plot-assistant/{{ report_data[\'folder_id\'] }}"' not in source

LOGIN_PARTIAL = 'shared/partials/pro_login_button.html'
LOGIN_MARKER = '<div id="cr-pro-login-btn">STUB</div>'


def _layout_env(partials):
    return Environment(
        loader=ChoiceLoader([FileSystemLoader(str(TEMPLATES)), DictLoader(partials)]),
        autoescape=True,
    )


def _layout_context(is_web, is_default_user=False, authenticated=True):
    return {
        'is_web': is_web,
        'is_default_user': is_default_user,
        'current_user': SimpleNamespace(is_authenticated=authenticated, email='a@b.c', role='User'),
        'config': {'BANNER_TEXT': 'CRISPResso'},
        'get_flashed_messages': lambda **kwargs: [],
        'url_for': lambda endpoint: '/' + endpoint,
    }


@pytest.mark.parametrize('is_default_user, authenticated', [(True, True), (False, True), (False, False)])
def test_web_layout_includes_login_partial_once(is_default_user, authenticated):
    html = _layout_env({LOGIN_PARTIAL: LOGIN_MARKER}).get_template('layout.html').render(
        **_layout_context(True, is_default_user, authenticated))
    assert html.count('cr-pro-login-btn') == 1


def test_web_layout_renders_without_login_partial():
    html = _layout_env({}).get_template('layout.html').render(**_layout_context(True, True))
    assert 'cr-pro-login-btn' not in html
    assert '</body>' in html


def test_cli_layout_never_includes_login_partial():
    with_partial = _layout_env({LOGIN_PARTIAL: LOGIN_MARKER}).get_template('layout.html').render(**_layout_context(False))
    without_partial = _layout_env({}).get_template('layout.html').render(**_layout_context(False))
    assert 'cr-pro-login-btn' not in with_partial
    assert with_partial == without_partial


def test_login_partial_sees_layout_context():
    """C2Web's partial reads current_user and is_default_user, so the include must pass the render context."""
    partial = '{% if current_user.email == "a@b.c" and is_default_user %}CTX-OK{% endif %}'
    html = _layout_env({LOGIN_PARTIAL: partial}).get_template('layout.html').render(**_layout_context(True, True))
    assert 'CTX-OK' in html
