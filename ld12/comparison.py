# Built-in modules #

# Internal modules #
from ld12 import families

# First party modules #
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class Comparison(object):
    """Takes care of comparing the trees created by the various clusters with
    the master ribosomal protein trees.

    For every tree (from clusters) compare with the ribosomal tree: How well does the structure correspond. How often do we find the same split between the families.When we don't fin the same split, which families are miss-spiting.

    1) Check that on the cluster trees the families are grouping together, either kick the ones that don't out, or report how many times a certain split is not conserved (uncollapsible)

    2) The trees that collapse nicely compare it with master tree setting the root at V and checking that the same branching patterns exists.

    3) If there is a popular alternative branching pattern, draw it in the supplementary and add as numbers how many clusters had this exact alternative branching."""

    def __init__(self, analysis):
        # Save attributes #
        self.analysis = analysis
        # The reference tree #
        self.analysis.ribosomal.master.tree_ete

    @property
    def ref_tree(self, analysis): return self.analysis.ribosomal.master.tree_dp

    @property_cached
    def collapsible(self):
        """Clusters that have strict coherence within all the seven families
        (i.e. all genes from every given family are monophyletic)"""
        # We are going to maintain two lists #
        collapsible = []
        uncollapsible = []
        # Check every one of the good clusters #
        for c in self.analysis.best_clusters:
            for f in families.values():
                if len([g for g in c if g.genome.family == f]) == 1: continue
                if c.tree_ete.check_monophyly(values=[f.name], target_attr="family")[0]: continue
                uncollapsible += c
                break
            else: collapsible += c
        # We are interested only in the collapsible #
        return collapsible

    @property_cached
    def matching(self, analysis):
        """Clusters that are collapsible and that have the same topology
        as the ribosomal master tree given that we collapse each family
        into one leaf."""
        # We are going to maintain two lists #
        matching = []
        mismatching = []
        # Check every one of the collapsible clusters #
        for c in self.collapsible:
            # Make a copy #
            tree = c.tree.copy(method='deepcopy')
            # Collapse families into one leaf #
            for f in families.values():
                fam_node = tree.get_common_ancestor(tree.search_nodes(family=f.name))
                for leaf in fam_node: assert leaf.family == f.name
                for child in fam_node.children: fam_node.remove_child(child)
                fam_node.name = f.name
            # Compare with rib master #
            rf, max_rf, common_leaves, parts_t1, parts_t2 = tree.robinson_foulds()
            if max_rf == 0: matching += c
            else:           mismatching += c
        # We are interested only in the matching #
        return matching