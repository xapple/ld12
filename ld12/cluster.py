# Built-in modules #

# First party modules #
from fasta import FASTA, AlignedFASTA
from plumbing.autopaths import FilePath

# Third party modules #
import ete2
from Bio import Phylo

###############################################################################
class Cluster(object):
    """A set of genes which are related in some way. For instance, all genes
    that clustered together when performing the MCL analysis.
    (all variations of a some single copy gene class)."""

    def __repr__(self): return '<%s object number %i>' % (self.__class__.__name__, self.num)

    def __init__(self, num, ids, analysis, name=None):
        # Basic params #
        self.num = num
        self.ids = ids
        self.analysis = analysis
        self.name = "cluster_%i" % num if name is None else name
        self.path = self.analysis.p.clusters_dir + self.name + '.fasta'
        # Pick only one gene per genome #
        self.genomes = [g for g in self.analysis.genomes if g.ids&self.ids]
        self.filtered_ids = [set(g.ids&self.ids).pop() for g in self.genomes]

    @property
    def counts(self):
        return self.analysis.count_table.loc[self.name]

    @property
    def score(self):
        """Given the genome counts, what is the single-copy likelihood score"""
        score = 0
        for name, count in self.counts.iteritems():
            partial = [g for g in self.analysis.genomes if g.prefix == name][0].partial
            if count == 0:   score +=  -5 if partial else -20
            elif count == 1: score +=  10 if partial else  10
            elif count == 2: score += -35 if partial else -30
            elif count == 3: score += -45 if partial else -40
            else: score += -100
        return score

    @property
    def sequences(self):
        for id_num in self.ids: yield self.analysis.blast_db[id_num]

    @property
    def fasta(self):
        """The fasta file containing the sequences of this cluster"""
        fasta = FASTA(self.path)
        if not fasta:
            with fasta as handle:
                for genome in self.genomes:
                    seq_id = set(genome.ids & self.ids).pop()
                    seq = self.analysis.blast_db[seq_id]
                    handle.add_str(str(seq.seq), name=genome.prefix)
        return fasta

    @property
    def alignment(self):
        """The fasta file aligned with muscle and filtered with gblocks"""
        alignment = AlignedFASTA(self.fasta.prefix_path + '.aln')
        if not alignment.exists:
            muscle = AlignedFASTA(self.fasta.prefix_path + '.muscle')
            self.fasta.align(muscle)
            print muscle.gblocks(alignment, nucleotide=self.analysis.sequence_type=='nucleotide')
        return alignment

    @property
    def tree(self):
        """The tree built with raxml"""
        tree = FilePath(self.alignment.prefix_path + '.tree')
        if not tree.exists:
            print "Building tree for cluster '%s'..." % self.name
            self.alignment.build_tree(tree, nucleotide=self.analysis.sequence_type=='nucleotide')
        return tree

    @property
    def phylogeny(self):
        """We can parse it with biopython"""
        return Phylo.read(self.tree.path, 'newick')

    def draw_tree(self):
        """120 pixels per branch length unit"""
        t = ete2.Tree(self.tree.contents)
        t.render(self.tree.replace_extension('pdf'))
