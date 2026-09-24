import json
import shlex
import sys

import pandas as pd
import pytest

from CRISPResso2 import CRISPRessoBatchCORE


def test_should_plot_large_plots():
    num_rows = 60
    c2pro_installed = False
    use_matplotlib = False
    large_plot_cutoff = 300
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_c2pro_installed_use_matplotlib_small():
    num_rows = 60
    c2pro_installed = True
    use_matplotlib = True
    large_plot_cutoff = 300
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_c2pro_installed():
    num_rows = 6000
    c2pro_installed = True
    use_matplotlib = False
    large_plot_cutoff = 300
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_c2pro_installed_use_matplotlib_large():
    num_rows = 6000
    c2pro_installed = True
    use_matplotlib = True
    large_plot_cutoff = 300
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_c2pro_not_installed_use_matplotlib():
    num_rows = 6000
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


# =============================================================================
# Additional edge case tests
# =============================================================================


def test_should_plot_large_plots_zero_rows():
    """Test with zero rows - should always plot."""
    num_rows = 0
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    # 0/6 = 0 < 300, so should plot
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_exact_cutoff_samples():
    """Test at exact cutoff boundary (in samples, not rows)."""
    # The function divides num_rows by 6 to get samples
    # 1800 rows / 6 = 300 samples, at cutoff boundary
    num_rows = 1800
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    # At exact cutoff, (1800/6 = 300) is NOT less than 300, so should not plot
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_one_below_cutoff_samples():
    """Test one sample below cutoff."""
    # 1794 rows / 6 = 299 samples, below cutoff
    num_rows = 1794
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_one_above_cutoff_samples():
    """Test one sample above cutoff."""
    # 1806 rows / 6 = 301 samples, above cutoff
    num_rows = 1806
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    # Above cutoff, should not plot (301 >= 300)
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_very_large_cutoff():
    """Test with very large cutoff."""
    num_rows = 10000
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 100000
    # 10000/6 = 1666.67 < 100000, should plot
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_small_cutoff():
    """Test with cutoff of 1 and enough rows to exceed it."""
    # 12 rows / 6 = 2 samples >= 1, should not plot
    num_rows = 12
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 1
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_c2pro_overrides_matplotlib():
    """Test that c2pro being installed overrides matplotlib restriction."""
    num_rows = 6000
    c2pro_installed = True
    use_matplotlib = False  # Not using matplotlib
    large_plot_cutoff = 300
    # Should plot because c2pro is installed and not using matplotlib
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_small_with_matplotlib():
    """Test small dataset with matplotlib."""
    num_rows = 10
    c2pro_installed = False
    use_matplotlib = True
    large_plot_cutoff = 300
    # 10/6 = 1.67 < 300, should plot
    assert CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


def test_should_plot_large_plots_not_using_matplotlib_large_no_c2pro():
    """Test large dataset not using matplotlib but c2pro not installed."""
    num_rows = 6000
    c2pro_installed = False
    use_matplotlib = False
    large_plot_cutoff = 300
    # Without c2pro: (not use_matplotlib and c2pro_installed) is False
    # And (6000/6 = 1000) >= 300, so should NOT plot
    assert not CRISPRessoBatchCORE.should_plot_large_plots(num_rows, c2pro_installed, use_matplotlib, large_plot_cutoff)


@pytest.mark.parametrize('source', ['none', 'config_file', 'config_json'])
def test_batch_config_inheritance_and_row_override(tmp_path, monkeypatch, source):
    """Build real Batch commands, then parse each child's config as Core does."""
    shared = CRISPRessoBatchCORE.CRISPRessoShared
    if source == 'config_json' and not shared.is_C2Pro_installed():
        pytest.skip('Inline config arguments require CRISPRessoPro')
    configs = [
        {'colors': {'A': '#123456'}, 'figures': [{'section_name': 'Alignment statistics', 'content': []}]},
        {'colors': {'A': '#654321'}, 'figures': [{'section_name': "O'Brien's alignment", 'content': []}]},
    ]
    values = [json.dumps(config) for config in configs]
    if source == 'config_file':
        paths = [tmp_path / "global config's.json", tmp_path / 'row config.json']
        for path, value in zip(paths, values):
            path.write_text(value)
        values = [str(path) for path in paths]
    fastq = tmp_path / 'reads.fastq'
    fastq.touch()
    rows = [
        {'name': name, 'fastq_r1': str(fastq), 'amplicon_seq': 'ACGT' * 25}
        for name in ['inherited', 'overridden']
    ]
    argv = ['CRISPRessoBatch', '-bs', str(tmp_path / 'batch.tsv'), '-bo', str(tmp_path)]
    if source != 'none':
        rows[0][source] = None
        rows[1][source] = values[1]
        argv += ['--' + source, values[0]]
    pd.DataFrame(rows).to_csv(tmp_path / 'batch.tsv', sep='\t', index=False)
    monkeypatch.setattr(sys, 'argv', argv)
    commands = []

    def capture_commands(cmds, *args, **kwargs):
        commands.extend(cmds)
        # Stop before analysis, without being swallowed by main's error handler.
        raise SystemExit(42)

    monkeypatch.setattr(CRISPRessoBatchCORE.CRISPRessoMultiProcessing, 'run_crispresso_cmds', capture_commands)
    old_handlers = list(CRISPRessoBatchCORE.logger.handlers)
    try:
        with pytest.raises(SystemExit) as exit_info:
            CRISPRessoBatchCORE.main()
        assert exit_info.value.code == 42
    finally:
        for handler in list(CRISPRessoBatchCORE.logger.handlers):
            if handler not in old_handlers:
                CRISPRessoBatchCORE.logger.removeHandler(handler)
                handler.close()
    assert len(commands) == 2
    for index, command in enumerate(commands):
        child_args = shared.getCRISPRessoArgParser('Core').parse_args(shlex.split(command)[1:])
        if source == 'none':
            assert child_args.config_file in (None, 'None')
            assert getattr(child_args, 'config_json', None) in (None, 'None')
        else:
            assert getattr(child_args, source) == values[index]
            if shared.is_C2Pro_installed():
                loaded = shared.check_custom_config(child_args)
                assert loaded['colors']['A'] == configs[index]['colors']['A']
                assert loaded['figures'] == configs[index]['figures']
