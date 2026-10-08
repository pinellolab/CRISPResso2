"""Tests that CLI startup does not import plotting libraries."""

import subprocess
import sys


def _run_isolated(code):
    result = subprocess.run(
        [sys.executable, '-c', code],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(
            f'isolated python failed (rc={result.returncode})\n'
            f'stdout:\n{result.stdout}\nstderr:\n{result.stderr}'
        )
    return result


def test_package_import_does_not_load_matplotlib():
    """import CRISPResso2 must not pull in matplotlib or seaborn."""
    _run_isolated(
        'import sys\n'
        'import CRISPResso2\n'
        'heavy = [m for m in ("matplotlib", "seaborn", "matplotlib.pyplot") if m in sys.modules]\n'
        'assert heavy == [], heavy\n'
    )


def test_core_import_does_not_load_matplotlib():
    """Importing CRISPRessoCORE must not pull in matplotlib or seaborn."""
    _run_isolated(
        'import sys\n'
        'from CRISPResso2.CRISPRessoCORE import main\n'
        'heavy = [m for m in ("matplotlib", "seaborn") if m in sys.modules]\n'
        'assert heavy == [], heavy\n'
        'assert callable(main)\n'
    )


def test_legacy_crispressoplot_import_path():
    """from CRISPResso2.CRISPRessoPlot import X still works for CRISPRessoPro."""
    _run_isolated(
        'from CRISPResso2.CRISPRessoPlot import get_nuc_color\n'
        'color = get_nuc_color("A", 1.0)\n'
        'assert abs(color[0] - 127 / 255.0) < 1e-9\n'
    )


def test_legacy_package_attribute_import():
    """from CRISPResso2 import CRISPRessoPlot still works."""
    _run_isolated(
        'from CRISPResso2 import CRISPRessoPlot\n'
        'assert callable(CRISPRessoPlot.get_nuc_color)\n'
    )


def test_cli_version_does_not_load_core_or_matplotlib():
    """CRISPResso --version must not import CORE or matplotlib."""
    result = _run_isolated(
        'import sys\n'
        'sys.argv = ["CRISPResso", "--version"]\n'
        'from CRISPResso2.cli import main\n'
        'rc = main()\n'
        'assert rc == 0\n'
        'heavy = [m for m in ("matplotlib", "seaborn", "pandas", "CRISPResso2.CRISPRessoCORE") if m in sys.modules]\n'
        'assert heavy == [], heavy\n'
    )
    assert 'CRISPResso 2.' in result.stdout


def test_cli_version_output():
    """--version prints the program name and version on stdout."""
    from CRISPResso2.startup import __version__

    result = subprocess.run(
        [sys.executable, '-c',
         'import sys; sys.argv = ["CRISPResso", "--version"]; '
         'from CRISPResso2.cli import main; raise SystemExit(main())'],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert result.stdout.strip() == f'CRISPResso {__version__}'
    assert result.stderr == ''


def test_cli_help_does_not_load_core_or_pandas():
    """CRISPResso --help must not import CORE, pandas, or matplotlib."""
    result = subprocess.run(
        [
            sys.executable,
            '-c',
            'import sys\n'
            'sys.argv = ["CRISPResso", "--help"]\n'
            'from CRISPResso2.cli import main\n'
            'try:\n'
            '    main()\n'
            'except SystemExit as exc:\n'
            '    assert exc.code in (0, None)\n'
            'heavy = [m for m in ("matplotlib", "seaborn", "pandas", "CRISPResso2.CRISPRessoCORE") if m in sys.modules]\n'
            'assert heavy == [], heavy\n',
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise AssertionError(result.stdout + result.stderr)
    assert '~~~CRISPResso 2~~~' in result.stdout
    assert 'usage:' in result.stdout.lower()
    first_content = next(line for line in result.stdout.splitlines() if line.strip())
    assert 'CRISPResso' in first_content
