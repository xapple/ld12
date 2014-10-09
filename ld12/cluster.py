# Built-in modules #
from collections import defaultdict

# Internal modules #
from ld12 import genes

# First party modules #
from fasta import FASTA, AlignedFASTA
from plumbing.autopaths import AutoPaths, FilePath

# Third party modules #
import ete2
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
    """

    def __repr__(self): return '<%s object number %i>' % (self.__class__.__name__, self.num)

    def __init__(self, num, line, analysis, name=None):
        # Basic params #
        self.num = num
        self.analysis = analysis
        # All the genes #
        self.genes = [genes[x] for x in frozenset(line.split())]
        # Optional #
        self.name = "cluster_%i" % num if name is None else name
        # Paths #
        self.base_dir = self.analysis.p.clusters_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Genomes #
        self.genomes = defaultdict(int)
        for gene in self.genes: self.genomes[gene.genome] += 1
        # Only one gene per genome #
        self.filtered_genes = [(set(G.genes.values()) & set(self.genes)).pop() for G in self.genomes]

    @property
    def counts(self):
        return self.analysis.count_table.loc[self.name]

    @property
    def score(self):
        """Given the genome counts, what is the single-copy custom likelihood score"""
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
            for gene in self.filtered_genes: fasta.add_str(str(gene), name=gene.long_name)
            fasta.close()
        return fasta

    @property
    def alignment(self):
        """The fasta file aligned with muscle and filtered with gblocks"""
        muscle    = AlignedFASTA(self.p.muscle)
        alignment = AlignedFASTA(self.p.aln)
        if not alignment:
            self.fasta.align(muscle)
            muscle.gblocks(self.p.aln, seq_type = self.analysis.seq_type)
        return alignment

    @property
    def tree(self):
        """The tree built with raxml"""
        tree = FilePath(self.alignment.prefix_path + '.tree')
        if not tree.exists:
            print "Building tree for cluster '%s'..." % self.name
            self.alignment.build_tree(tree,
                                      seq_type = self.analysis.seq_type,
                                      num_threads = self.analysis.num_threads)
        return tree

    @property
    def phylogeny(self):
        """We can parse it with biopython"""
        return Phylo.read(self.tree.path, 'newick')

    def draw_tree(self):
        """120 pixels per branch length unit"""
        t = ete2.Tree(self.tree.contents)
        t.render(self.tree.replace_extension('pdf'))
