# Built-in modules #

# Internal modules #
from ld12 import genomes

# First party modules #
from seqsearch.blast import BLASTdb
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from seqsearch.databases.refseq import refseq_bact_prot_nr, refseq_arch_prot_nr

# Third party modules #
from shell_command import shell_output

###############################################################################
class RefSeqProkPlusMarine(object):
    """A special database combining the refseq non redundant protein archives
    for both bacteria and archaea, with, as an extra, all the marine genes added.
    In addition, we will also add the freshwater genomes that appear to have been
    already added to RefSeq."""

    all_paths = """
    /all_genes.fasta
    /all_genes.fasta.00.pin
    /log.txt
    /out.txt
    """

    def __init__(self, base_dir, duplications):
        # Base attributes #
        self.duplications = duplications
        self.timer = duplications.timer
        # Directories #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra parameters #
        self.refseq_bact_orig = list(refseq_bact_prot_nr.raw_files)
        self.refseq_arch_orig = list(refseq_arch_prot_nr.raw_files)
        self.missing_marine   = [g for g in genomes.values() if g.marine and not g.in_refseq_bact]
        self.missing_fresh    = [g for g in genomes.values() if g.fresh and not g.in_refseq_bact]

    @property_cached
    def blast_db(self):
        """A blastable database of refseq + all marine organism genes"""
        blast_db = BLASTdb(self.p.genes, 'prot')
        if not self.p.genes.exists:
            # We are going to cat a whole of files together #
            print "Regrouping all fasta files together..."
            all_genes = self.refseq_bact_orig   \
                      + self.refseq_arch_orig \
                      + self.missing_marine   \
                      + self.missing_fresh
            shell_output("zcat %s > %s" % (' '.join(all_genes), self.p.genes))
            self.timer.print_elapsed()
            # Check that all files ended with a newline #
            print "Checking that sequence counts match..."
            assert len(blast_db) == sum(map(len,self.refseq_bact_orig))  \
                                  + sum(map(len,self.refseq_arch_orig)) \
                                  + sum(map(len,self.missing_marine))   \
                                  + sum(map(len,self.missing_fresh))
            self.timer.print_elapsed()
        if not self.p.pin.exists:
            # Call make DB #
            print "Building a BLAST database..."
            blast_db.makeblastdb(logfile=self.p.log, out=self.p.out)
            self.timer.print_elapsed()
        return blast_db
