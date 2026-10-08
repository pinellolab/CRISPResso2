# Backward-compatible re-exports: CRISPRessoPlot and upsetplot moved to
# CRISPResso2.plots but can still be imported from the top-level package.
# We can remove this once we update CRISPRessoPro's importing.
#
# These aliases are lazy so `import CRISPResso2` (and therefore every CLI
# entry point) does not pay for matplotlib/seaborn until plots are used.
import importlib
import sys
import types


class _LazyModule(types.ModuleType):
    """Stand-in for a submodule that is imported on first attribute access."""

    def __init__(self, name, target):
        super().__init__(name)
        self.__dict__['_target'] = target
        self.__dict__['_loaded'] = None

    def _load(self):
        loaded = self.__dict__['_loaded']
        if loaded is None:
            loaded = importlib.import_module(self.__dict__['_target'])
            self.__dict__['_loaded'] = loaded
            sys.modules[self.__name__] = loaded
            parent_name, _, attr = self.__name__.rpartition('.')
            parent = sys.modules.get(parent_name)
            if parent is not None and getattr(parent, attr, None) is self:
                setattr(parent, attr, loaded)
        return loaded

    def __getattr__(self, name):
        return getattr(self._load(), name)

    def __dir__(self):
        return dir(self._load())


def _install_lazy_alias(alias, target):
    """Register `alias` in sys.modules as a lazy wrapper around `target`."""
    lazy = _LazyModule(alias, target)
    sys.modules[alias] = lazy
    return lazy


# Register under old module paths so that `from CRISPResso2.CRISPRessoPlot import X`
# and `from CRISPResso2.upsetplot import X` still work for downstream consumers
# (e.g., CRISPRessoPro). Remove once all consumers are updated.
CRISPRessoPlot = _install_lazy_alias(
    'CRISPResso2.CRISPRessoPlot',
    'CRISPResso2.plots.CRISPRessoPlot',
)
upsetplot = _install_lazy_alias(
    'CRISPResso2.upsetplot',
    'CRISPResso2.plots.upsetplot',
)
