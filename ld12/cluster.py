# Built-in modules #
from collections import defaultdict

# Internal modules #
from ld12 import genes

# First party modules #
from fasta import FASTA, AlignedFASTA
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache import property_cached

# Third party modules #
import ete2, dendropy
from Bio import Phylo

###############################################################################
class Cluster(object):
    """A set of genes which are related in some way. For instance, all genes
    that clustered together when performing the MCL analysis.
    (e.g. all variants of a some single copy gene class)."""

    all_paths = """
    /filtered_genes.fasta
    /filtered_genes.muscle
    /filtered_genes.aln
    /tree/
    /tree/RAxML_bestTree.tree
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __len__(self): return len(self.genes)
    def __iter__(self): return iter(self.genes)

    def __init__(self, num, line, analysis, name=None):
        # Basic params #
        self.num = num
        self.analysis = analysis
        # All the genes #
        self.genes = [genes[x] for x in frozenset(line.split())]
        # Optional #
        self.name = "cluster_%s" % num if name is None else name
        # Paths #
        self.base_dir = self.analysis.p.clusters_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Genomes #
        self.genomes = defaultdict(int)
        for gene in self.genes: self.genomes[gene.genome] += 1
        # Only one gene per genome #
        self.filtered_genes = [(set(G.genes.values()) & set(self.genes)).pop() for G in self.genomes]

    @property
    def families(self):
        """Which families are represented by this cluster"""
        return set([g.genome.family for g in self.filtered_genes])

    @property
    def counts(self):
        """Before filtering, how many genes in each genomes for this cluster"""
        return self.analysis.count_table[self.name]

    @property
    def score(self):
        """Given the genome counts, what is our single-copy custom likelihood score"""
        score = 0
        for genome, count in self.genomes.items():
            if   count == 0: score +=  -5 if genome.partial else -20
            elif count == 1: score +=  10 if genome.partial else  10
            elif count == 2: score += -35 if genome.partial else -30
            elif count == 3: score += -45 if genome.partial else -40
            else:            score += -100
        return score

    @property
    def fasta(self):
        """The fasta file containing the filtered genes of this cluster
        The names now will correspond to long descriptive names"""
        fasta = FASTA(self.p.fasta)
        if not fasta:
            fasta.create()
            for gene in self.filtered_genes: fasta.add_str(str(gene), name=gene.name)
            fasta.close()
        return fasta

    @property
    def alignment(self):
        """The fasta file aligned with muscle and filtered with gblocks"""
        muscle    = AlignedFASTA(self.p.muscle)
        alignment = AlignedFASTA(self.p.aln)
        if not alignment:
            self.fasta.align(muscle)
            muscle.gblocks(self.p.aln, seq_type=self.analysis.seq_type)
        return alignment

    @property
    def tree(self):
        """The path to the tree built with raxml"""
        tree = FilePath(self.p.tree_dir + 'RAxML_bestTree.tree')
        if not tree.exists:
            self.alignment.build_tree(new_path    = self.p.tree_dir,
                                      seq_type    = self.analysis.seq_type,
                                      num_threads = self.analysis.num_threads,
                                      free_cores  = 0,
                                      keep_dir    = True)
        return tree

    @property_cached
    def tree_dp(self):
        """The tree as an object in python memory from dendropy
        We can modify the leaves to link to their Gene object"""
        return dendropy.Tree.get_from_path(self.tree, 'newick')

    @property_cached
    def tree_biop(self):
        """The tree as an object in python memory from biopython"""
        return Phylo.read(self.tree, 'newick')

    @property_cached
    def tree_ete(self):
        """The tree as an object in python memory from ETE2
        We can add attributes to the leaves useful for the comparisons
        that we perform later on."""
        # Load it #
        tree = ete2.Tree(self.tree)
        # Add info to leaves #
        for leaf in tree:
            gene = genes[leaf.name]
            leaf.add_features(gene=gene)
            leaf.add_features(family=gene.genome.family.name)
        # Root it #
        five = tree.search_nodes(family='v')
        assert len(five) == 1
        tree.set_outgroup(five[0])
        tree.ladderize()
        # Return results #
        return tree

    def print_tree(self):
        """Render it as ASCII on your terminal"""
        print self.tree_dp.as_ascii_plot()

    def draw_tree(self):
        """120 pixels per branch length unit"""
        self.tree_ete.render(self.tree.replace_extension('pdf'))
