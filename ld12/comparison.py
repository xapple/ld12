# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import dendropy

###############################################################################
class Comparison(object):
    """Takes care of comparing the trees created by the various clusters with
    the master ribosomal protein trees.

    For every tree (from clusters) compare with the ribosomal tree: How well does the structure correspond. How often do we find the same split between the families.When we don't fin the same split, which families are miss-spiting.

    1) Check that on the cluster trees the families are grouping together, either kick the ones that don't out, or report how many times a certain split is not conserved (uncollapsible)

    2) The trees that collapse nicely compare it with master tree setting the root at V and checking that the same branching patterns exists.

    3) If there is a popular alternative branching pattern, draw it in the supplementary and add as numbers how many clusters had this exact alternative branching."""

    def __init__(self, analysis):
        self.analysis = analysis

    @property
    def ref_tree(self, analysis): return self.analysis.ribosomal.master.tree_dp

    @property_cached
    def collapsable(self, analysis):
        """Clusters that have strict coherence within all the seven families
        (i.e. all genes from every given family are monophyletic)"""
        for c in self.analysis.best_clusters:
            pass

    @property_cached
    def matching(self, analysis):
        """Clusters that have the"""
        pass

    @property_cached
    def mismatching(self, analysis):
        pass
