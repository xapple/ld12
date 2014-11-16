# Built-in modules #

# Internal modules #
from ld12 import families

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import pandas

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

    @property
    def ref_tree(self, analysis): return self.analysis.ribosomal.master.tree_dp

    @property_cached
    def collapsing(self):
        """Clusters that have strict coherence within all the seven families
        (i.e. all genes from every given family are monophyletic) are
        considered collapsible"""
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
        # Return values #
        return collapsible, uncollapsible

    @property
    def collapsible(self): return self.collapsing[0]
    @property
    def uncollapsible(self): return self.collapsing[1]

    @property_cached
    def matches(self, analysis):
        """Taking only clusters that are collapsible, we can ask: do they
        match the same topology as the ribosomal master tree ? This can be
        done since we can collapse each family into one leaf."""
        # We are going to maintain two lists #
        matching = []
        mismatching = []
        # Check every one of the collapsible clusters #
        for c in self.collapsible:
            # Make a copy #
            tree = c.tree_ete.copy(method='deepcopy')
            # Collapse families into one leaf #
            for f in families.values():
                nodes = tree.search_nodes(family=f.name)
                if len(nodes) == 1: ancestor = nodes[0]
                else:               ancestor = tree.get_common_ancestor(nodes)
                for leaf in ancestor: assert leaf.family == f.name
                for node in ancestor.iter_descendants(): node.delete(prevent_nondicotomic=False)
                ancestor.name = f.name
            # Compare with rib master reference #
            ref = self.analysis.ribosomal.master.tree_ete
            rf, max_rf, common_leaves, parts_t1, parts_t2 = tree.robinson_foulds(ref)
            print c, rf, max_rf, common_leaves, parts_t1, parts_t2
            if rf == 0: matching += c
            else:       mismatching += c
        # Return values #
        return matching, mismatching

    @property
    def matching(self): return self.matches[0]
    @property
    def mismatching(self): return self.matches[1]

    #-------------------------------------------------------------------------#
    def uncollapsible_stats(self):
        """For the trees that are uncollapsible, where do they not correspond?"""
        result = {}
        for c in self.uncollapsible:
            result[c.name] = {}
            for f in families.values():
                if len([g for g in c if g.genome.family == f]) == 1:
                    result[f.name] = 'single'
                elif c.tree_ete.check_monophyly(values=[f.name], target_attr="family")[0]:
                    result[f.name] = 'mono'
                else:
                    tree = c.tree_ete
                    nodes = tree.search_nodes(family=f.name)
                    ancestor = tree.get_common_ancestor(nodes)
                    intruders = default_dict(int)
                    for leaf in ancestor: intruders[leaf.family] += 1
                    intruders.pop(f.name)
                    result[f.name] = dict(intruders)
        result = pandas.DataFrame(result)
        return result

    def save_uncollapsible_stats(self):
        """Save the dataframe above in a CSV file"""
        self.uncollapsible_stats.to_csv(str(self.analysis.p.uncollapsible), sep='\t', encoding='utf-8')

    def mismatching_stats(self):
        """For the trees that mismatch the reference, where do they not
        correspond?"""
        result = {}
        for c in self.uncollapsible:
            pass

    def save_mismatching_stats(self):
        """Save the dataframe above in a CSV file"""
        self.mismatching_stats.to_csv(str(self.analysis.p.mismatching), sep='\t', encoding='utf-8')