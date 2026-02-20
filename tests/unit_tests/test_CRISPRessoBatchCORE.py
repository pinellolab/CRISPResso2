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
