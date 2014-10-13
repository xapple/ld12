# Built-in modules #
from collections import defaultdict

# Internal modules #
from ld12 import genes, families
from ld12.cluster import Cluster
from ld12.master import Master

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import pandas
from tqdm import tqdm

###############################################################################
class Ribosomal(object):
    """Takes care of generating the extra ribosomal protein clusters and trees."""

    def __repr__(self): return '<%s object with %i groups>' % \
        (self.__class__.__name__, len(self.groups))

    def __init__(self, analysis):
        self.analysis = analysis
        self.genes = [g for g in genes.values() if g.ribo_group]
        self.groups = set([g.ribo_group for g in self.genes])
        self.master = Master(self)

    @property_cached
    def clusters(self):
        """One cluster for every ribo type."""
        clusters = []
        for i, group in enumerate(self.groups):
            current_genes = [g for g in self.genes if g.ribo_group==group]
            current_genes = ' '.join(g.name for g in current_genes)
            name = "cluster_r%i" % i
            cluster = Cluster(i, current_genes, self.analysis, name)
            cluster.ribo_group = group
            clusters.append(cluster)
        clusters = sorted(clusters, key=lambda x: x.score, reverse=True)
        return clusters

    def make_trees(self):
        """If you access a tree it will be built, but as it takes time,
        let's all do them together now in a non-lazy way."""
        for i, c in enumerate(self.clusters):
            print "* Building tree for ribosomal cluster '%s'..." % c.name
            print "* Cluster %i out of %i of group '%s' with %i genes (%i filtered) and a score of %i" % \
                  (i+1, len(self.clusters)+1, len(c.genes), c.genes[0].ribo_group, len(c.filtered_genes), c.score)
            print c.tree
            self.analysis.timer.print_elapsed()

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

    #-------------------------------------------------------------------------#
    @property_cached
    def family_clusters_table(self):
        """A table evaluating how many ribo_type representative there are per family."""
        result = defaultdict(lambda: defaultdict(int))
        for group in self.groups:
            current_genes = [g for g in self.genes if g.ribo_group==group]
            for gene in current_genes:
                result[group][gene.genome.family.name] += 1
        result = pandas.DataFrame(result, dtype=int)
        result = result.fillna(0)
        result = result.astype(int)
        result = result.reindex_axis(result.sum(axis=1).order(ascending=False).index)
        result = result.reindex_axis(result.sum().order(ascending=False).index, axis='columns')
        return result

    def save_family_table(self):
        self.family_clusters_table.to_csv(str(self.analysis.p.family_tsv), sep='\t', encoding='utf-8')

    #-------------------------------------------------------------------------#
    @property_cached
    def family_clusters(self):
        """A special set of clusters combining fewer ribo groups with only one representative per family."""
        clusters = []
        i = 0
        for c in self.clusters:
            # Check there is one representative per family #
            if len(c.families) != len(families): continue
            # Pick one per family #
            current_genes = [(set(f.genes) & set(c.genes)).pop() for f in families.values()]
            # Make cluster #
            current_genes = ' '.join(g.name for g in current_genes)
            name = "cluster_f%i" % i
            cluster = Cluster(i, current_genes, self.analysis, name)
            cluster.ribo_group = c.ribo_group
            clusters.append(cluster)
            i += 1
        return clusters

    def make_alignments(self):
        """If you access a alignment it will be built, but as it takes time,
        let's all do them together now in a non-lazy way."""
        for c in tqdm(self.family_clusters):
            assert c.alignment.ids == set(g.name for g in c.genes)