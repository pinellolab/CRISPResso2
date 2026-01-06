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


