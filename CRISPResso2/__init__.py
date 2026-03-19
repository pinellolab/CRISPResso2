# Backward-compatible re-exports: CRISPRessoPlot and upsetplot moved to
# CRISPResso2.plots but can still be imported from the top-level package.
# We can remove this once we update CRISPRessoPro's importing
import sys

from CRISPResso2.plots import CRISPRessoPlot
from CRISPResso2.plots import upsetplot

# Register under old module paths so that `from CRISPResso2.CRISPRessoPlot import X`
# and `from CRISPResso2.upsetplot import X` still work for downstream consumers
# (e.g., CRISPRessoPro). Remove once all consumers are updated.
sys.modules['CRISPResso2.CRISPRessoPlot'] = CRISPRessoPlot
sys.modules['CRISPResso2.upsetplot'] = upsetplot
