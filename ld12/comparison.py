# Futures #
from __future__ import division

# Built-in modules #
import string
from collections import defaultdict

# Internal modules #
from ld12 import families

# First party modules #
from plumbing.cache import property_cached
from plumbing.common import pad_with_whitespace, mirror_lines, concatenate_by_line

# Third party modules #
import pandas
from tqdm import tqdm

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

    def collapse_families(self, tree, fams):
        """Utility function: Collapse families into one leaf"""
        for f in fams:
            nodes = tree.search_nodes(family=f.name)
            if len(nodes) == 1: ancestor = nodes[0]
            else:               ancestor = tree.get_common_ancestor(nodes)
            for leaf in ancestor: assert leaf.family == f.name
            for node in ancestor.iter_descendants(): node.delete(prevent_nondicotomic=False)
            ancestor.name = f.name
            ancestor.family = f.name

    #-------------------------------------------------------------------------#
    #                             COLLAPSING                                  #
    #-------------------------------------------------------------------------#
    @property_cached
    def collapsing(self):
        """Clusters that have strict coherence within all the seven families
        (i.e. all genes from every given family are monophyletic) are
        considered collapsible"""
        # We are going to maintain two lists #
        collapsible = []
        uncollapsible = []
        # Check every one of the good clusters #
        print "Computing which clusters are collapsing..."
        for c in tqdm(self.analysis.best_clusters):
            if not c.p.bestTree.exists:
                print "Warning: cluster %s is missing a tree, skipping" % c
                continue
            for f in families.values():
                if len([g for g in c if g.genome.family == f]) == 1: continue
                if c.tree_ete.check_monophyly(values=[f.name], target_attr="family")[0]: continue
                uncollapsible += [c]
                break
            else: collapsible += [c]
        # Return values #
        return collapsible, uncollapsible

    @property
    def collapsible(self): return self.collapsing[0]
    @property
    def uncollapsible(self): return self.collapsing[1]

    @property
    def uncollapsible_stats(self):
        """For the trees that are not collapsible, where do they not correspond?"""
        result = {}
        for c in self.uncollapsible:
            result[c.name] = {}
            for f in families.values():
                if len([g for g in c if g.genome.family == f]) == 1:
                    result[c.name][f.name] = 'single'
                elif c.tree_ete.check_monophyly(values=[f.name], target_attr="family")[0]:
                    result[c.name][f.name] = 'mono'
                else:
                    tree = c.tree_ete
                    nodes = tree.search_nodes(family=f.name)
                    ancestor = tree.get_common_ancestor(nodes)
                    intruders = defaultdict(int)
                    for leaf in ancestor: intruders[leaf.family] += 1
                    intruders.pop(f.name)
                    result[c.name][f.name] = dict(intruders)
        # Make a dataframe #
        result = pandas.DataFrame(result)
        # Calculate a summary #
        summary = {}
        for f in families.values():
            ok = ((result.loc[f.name] == 'mono') | (result.loc[f.name] == 'single')).sum()
            summary[f.name] = ok / len(self.uncollapsible)
        # Put everything together #
        result['summary'] = pandas.Series(summary)
        result = result.transpose()
        return result

    def save_uncollapsible_stats(self):
        """Save the dataframe above in a CSV file"""
        self.uncollapsible_stats.to_csv(str(self.analysis.p.uncollapsible), sep='\t', encoding='utf-8')

    #-------------------------------------------------------------------------#
    #                             SPLIT 3a-3b                                 #
    #-------------------------------------------------------------------------#
    @property_cached
    def split_three_a_b(self):
        """Do the trees that are monophyletic for IIIa and IIIb match the reference tree
        just for that split."""
        # Let's maintain two lists #
        split_conserved = []
        split_broken    = []
        # Families #
        fams = [f for name,f in families.items() if name == 'IIIa' or name == 'IIIb']
        # Trees that are monophyletic for IIIa and IIIb #
        good_clusters = ((self.uncollapsible_stats['IIIa'] == 'mono') | \
                         (self.uncollapsible_stats['IIIa'] == 'single')) & \
                        ((self.uncollapsible_stats['IIIb'] == 'mono') | \
                         (self.uncollapsible_stats['IIIb'] == 'single'))
        good_clusters = [c for c in self.uncollapsible if good_clusters[c.name]]
        # These clusters also match the query #
        good_clusters += self.collapsible
        # Let's check if IIIa and IIIb are linked by a node #
        for cluster in good_clusters:
            tree = cluster.tree_ete.copy(method='deepcopy')
            self.collapse_families(tree, fams)
            a = tree.search_nodes(family='IIIa')
            b = tree.search_nodes(family='IIIb')
            assert len(a) == 1
            assert len(b) == 1
            # Do they share the same parent #
            if not (a[0].up is b[0].up):
                split_broken += [cluster]
                continue
            # Is the bootstrap value sufficient #
            if a[0].support >= 80:
                split_conserved += [cluster]
            else:
                split_broken    += [cluster]
        # Results #
        return split_conserved, split_broken

    def save_split_three_a_b(self):
        """Save the dataframe above in a CSV file"""
        content = 'Conserved: %i\nBroken: %i' % (len(self.split_three_a_b[0]), len(self.split_three_a_b[1]))
        self.analysis.p.three.write(content)

    #-------------------------------------------------------------------------#
    #                          SPLITS CONSERVED                               #
    #-------------------------------------------------------------------------#
    @property_cached
    def splits_conserved(self):
        """For every tree, are the splits of the master tree conserved ?
        Is the tree monophyletic for 3a and 3b? Is the tree monophyletic
        for 3a, 3b and 5 ? etc. all the way up"""
        # The result for every tree #
        result = {}
        groups = [('IIIa', 'IIIb', 'V', 'II', 'Ic', 'Ia', 'Ib')]
        # Main loop #
        for c in tqdm(self.analysis.best_clusters):
            if not c.p.bestTree.exists:
                print "Warning: cluster %s is missing a tree, skipping" % c
                continue
            for i, g in enumerate(groups):
                group = groups[0:i+1]
                conserved = c.tree_ete.check_monophyly(values=[group], target_attr="family")[0]
                result[c.name][group] = conserved
        # Make a dataframe #
        result = pandas.DataFrame(result)
        # Calculate a summary #
        summary = {}
        for i, g in enumerate(groups):
            group = groups[0:i+1]
            ok = (result.loc[group] == True)
            summary[group] = ok / len(self.analysis.best_clusters)
        # Put everything together #
        result['summary'] = pandas.Series(summary)
        result = result.transpose()
        return result

    def save_split_conserved(self):
        """Save the dataframe above in a CSV file"""
        self.splits_conserved.to_csv(str(self.analysis.p.conserved), sep='\t', encoding='utf-8')

    #-------------------------------------------------------------------------#
    #                          MATCHING WITH REF                              #
    #-------------------------------------------------------------------------#
    @property_cached
    def matches(self):
        """Taking only clusters that are collapsible, we can ask: do they
        match the same topology as the ribosomal master tree ? This can be
        done since we can collapse each family into one leaf."""
        # We are going to maintain two lists and one string #
        matching = []
        mismatching = []
        mismatching_stats = ""
        # Check every one of the collapsible clusters #
        print "Computing which clusters are matching..."
        for c in tqdm(self.collapsible):
            # Make a copy #
            tree = c.tree_ete.copy(method='deepcopy')
            # Collapse families into one leaf #
            self.collapse_families(tree, families.values())
            # Compare with rib master reference #
            ref = self.analysis.ribosomal.master.tree_ete
            rf, max_rf, common_leaves, parts_t1, parts_t2 = tree.robinson_foulds(ref)
            if rf == 0: matching += [c]
            else:
            # Let's collect some mismatching statistics #
                mismatching += [c]
                mismatching_stats += 'Tree from cluster %i is mismatching:\n' % c.num
                ref_string  = ref.get_ascii(show_internal=False)
                tree_string = tree.get_ascii(show_internal=False)
                ref_string  = pad_with_whitespace(ref_string)
                tree_string = pad_with_whitespace(tree_string)
                tree_string = mirror_lines(tree_string)
                tree_string = tree_string.translate(string.maketrans("/\\", "\\/"))
                mismatching_stats += concatenate_by_line(ref_string, tree_string) + '\n\n'
                mismatching_stats += "Robinson-Foulds metric: %i\n" % rf
                mismatching_stats += "Max RF: %i\n" % max_rf
                mismatching_stats += "-----------------------------------------------\n"
        # Return values #
        return matching, mismatching, mismatching_stats

    @property
    def matching(self): return self.matches[0]
    @property
    def mismatching(self): return self.matches[1]
    @property
    def mismatching_stats(self): return self.matches[2]

    def save_mismatching_stats(self):
        """Save the dataframe above in a CSV file"""
        self.analysis.p.mismatching.writelines(self.mismatching_stats)