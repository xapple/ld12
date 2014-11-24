# Built-in modules #

# Internal modules #
from ld12 import genomes

# First party modules #
from seqsearch.blast import BLASTdb
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta.databases import refseq_bact_prot_nr, refseq_arch_prot_nr

# Third party modules #
from shell_command import shell_output

###############################################################################
class RefSeqProkPlusMarine(object):
    """A special database combining the refseq non redundant protein archives
    for both bacteria and archaea, with, as an extra, all the marine genes added."""

    all_paths = """
    /all_genes.fasta
    /all_genes.fasta.nin
    /log.txt
    """

    def __init__(self, base_dir):
        # Base attributes #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra parameters #
        self.refseq_bact    = list(refseq_bact_prot_nr.raw_files)
        self.refseq_arch    = list(refseq_arch_prot_nr.raw_files)
        self.marine_genomes = [g for g in genomes.values() if g.marine]
        self.all_genes      = self.refseq_bact + self.refseq_arch + self.marine_genomes

    @property_cached
    def blast_db(self):
        """A blastable database of refseq + all marine organism genes"""
        blast_db = BLASTdb(self.p.genes, 'prot')
        if not self.p.genes.exists:
            # We are going to cat a whole of files together #
            shell_output("zcat %s > %s" % (' '.join(self.all_genes), self.p.fasta))
            # Check that all files ended with a newline #
            assert len(blast_db) == map(len,self.refseq_bact) + map(len,self.refseq_arch) + map(len,self.marine_genomes)
        if not self.p.nin.exists:
            # Call make DB #
            blast_db.makeblastdb(logfile=self.p.log)
        return blast_db
