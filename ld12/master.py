# Built-in modules #

# Internal modules #
from ld12.cluster import Cluster
from ld12 import families

# First party modules #
from fasta import AlignedFASTA
from plumbing.cache import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #
import ete2

###############################################################################
class Master(Cluster):
    """A master alignment and corresponding tree to be used to compare against
    other trees. A special cluster combining all the ribosomal family clusters in
        a concatenation."""

    def __init__(self, ribosomal):
        # Links
        self.ribosomal = ribosomal
        self.analysis = ribosomal.analysis
        # Basic params #
        self.num = -1
        self.name = "rib_master"
        # Paths #
        self.base_dir = self.analysis.p.clusters_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def alignment(self):
        alignment = AlignedFASTA(self.p.aln)
        if not alignment:
            alignment.create()
            for f in families.values():
                seq = ''
                for c in self.ribosomal.family_clusters:
                    gene = (set(f.genes) & set(c.genes)).pop()
                    assert gene.genome.family == f
                    assert gene.ribo_group == c.ribo_group
                    seq += str(c.alignment.sequences[gene.name].seq)
                alignment.add_str(seq, name=f.name)
            alignment.close()
        return alignment

    def make_tree(self):
        print "* Building tree for master ribosomal alignment"
        print self.tree
        self.analysis.timer.print_elapsed()

    @property_cached
    def tree_ete(self):
        """The tree as an object in python memory from ETE2
        We can add attributes to the leaves useful for the comparisons
        that we perform later on."""
        # Load it #
        tree = ete2.Tree(self.tree)
        # Root it #
        five = tree.search_nodes(name='v')
        assert len(five) == 1
        tree.set_outgroup(five[0])
        tree.ladderize()
        # Return results #
        return tree
