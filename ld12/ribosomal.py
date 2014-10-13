# Built-in modules #
from collections import defaultdict

# Internal modules #
from ld12 import genes
from ld12.cluster import Cluster

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class Ribosomal(object):
    """Takes care of generating the extra ribosomal protein clusters and trees."""

    def __repr__(self): return '<%s object with %i groups>' % \
        (self.__class__.__name__, len(self.groups))

    def __init__(self, analysis):
        self.analysis = analysis
        self.genes = [g for g in genes.values() if g.ribo_group]
        self.groups = set([g.ribo_group for g in self.genes])

    @property_cached
    def clusters(self):
        """One cluster for every ribo type."""
        clusters = []
        for i, group in enumerate(self.groups):
            current_genes = [g for g in self.genes if g.ribo_group==group]
            current_genes = ' '.join(g.name for g in current_genes)
            name = "cluster_r%i" % i
            clusters.append(Cluster(i, current_genes, self.analysis, name))
        return clusters

    @property_cached
    def family_cluster(self):
        """A special cluster with fewer ribo groups and only one representative per family."""
        clusters = []
        for i, group in enumerate(self.groups):
            current_genes = [g for g in self.genes if g.ribo_group==group]
            current_genes = ' '.join(g.name for g in current_genes)
            name = "r" + str(i)
            clusters.append(Cluster(name, current_genes, self.analysis))
        return clusters

    def make_trees(self):
        """If you access a tree it will be built, but as it takes time,
        let's all do them together now in a non-lazy way."""
        for i, c in enumerate(self.clusters):
            print "Building tree for ribosomal cluster '%s'..." % c.name
            print "Cluster %i out of %i with %i genes (%i filtered) and a score of %i" % \
                  (i+1, len(self.clusters)+1, len(c.genes), len(c.filtered_genes), c.score)
            print c.tree
            self.timer.print_elapsed()

    #-------------------------------------------------------------------------#
    @property_cached
    def frame(self):
        """A dataframe counting which ribo groups appear how many time in each genome."""
        result = defaultdict(lambda: defaultdict(int))
        for group in self.groups:
            current_genes = [g for g in self.genes if g.ribo_group==group]
            for gene in current_genes:
                result[group][gene.genome.name] += 1
        result = pandas.DataFrame(result, dtype=int)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(result.sum(axis=1).order(ascending=False).index)
        result = result.reindex_axis(result.sum().order(ascending=False).index, axis='columns')
        return result

    def save_ribo_table(self):
        self.frame.to_csv(str(self.analysis.p.ribo_tsv), sep='\t', encoding='utf-8')

    @property_cached
    def ref_tree(self):
        """A reference tree built with a number of ribosomal proteins (not 16S) to
        be used when comparing other trees to it."""
        pass