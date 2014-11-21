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
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def blast_db(self):
        """A blastable database of refseq + all marine organism genes"""
        blast_db = BLASTdb(self.p.genes, 'prot')
        if not self.p.genes.exists:
            # We are going to cat a whole of files together #
            refseq_bact  = list(refseq_bact_prot_nr.raw_files)
            refseq_arch  = list(refseq_arch_prot_nr.raw_files)
            marine_genomes = [g for g in genomes if g.marine]
            all_genes = refseq_bact + refseq_arch + marine_genomes
            shell_output("zcat %s > %s" % (' '.join(all_genes), self.p.fasta))
        if not self.p.nin.exists:
            # Call make DB #
            #blast_db.makeblastdb()
            pass
        return blast_db
