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
    /all_genes.fasta.00.pin
    /log.txt
    """

    def __init__(self, base_dir, duplications):
        # Base attributes #
        self.base_dir = base_dir
        self.duplications = duplications
        self.timer = duplications.timer
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
            print "Regrouping all fasta files together..."
            shell_output("zcat %s > %s" % (' '.join(self.all_genes), self.p.fasta))
            self.timer.print_elapsed()
            # Check that all files ended with a newline #
            print "Checking that sequence counts match..."
            assert len(blast_db) == map(len,self.refseq_bact) + map(len,self.refseq_arch) + map(len,self.marine_genomes)
            self.timer.print_elapsed()
        if not self.p.pin.exists:
            # Call make DB #
            print "Building a BLAST database..."
            blast_db.makeblastdb(logfile=self.p.log)
            self.timer.print_elapsed()
        return blast_db
