# Built-in modules #

# Internal modules #
from ld12 import genomes

# First party modules #
from seqsearch.blast import BLASTdb
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta.databases import refseq_bact_prot_nr, refseq_arch_prot_nr
from fasta import FASTA

# Third party modules #
from shell_command import shell_output

###############################################################################
class RefSeqProkPlusMarine(object):
    """A special database combining the refseq non redundant protein archives
    for both bacteria and archaea, with, as an extra, all the marine genes added.
    In addition, we will remove the freshwater genomes that appear to have been
    already added to RefSeq."""

    all_paths = """
    /modified_refseq_bact/
    /all_genes.fasta
    /all_genes.fasta.00.pin
    /log.txt
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
        self.marine_genomes   = [g for g in genomes.values() if g.marine]

    @property
    def refseq_bact_mod(self):
        """Copy the files from the original refseq bacteria database.
        Edit them to remove the things the genes don't want included
        and return a list of the new files."""
        # Message #
        print "Filtering refseq database by removing SCGC genes..."
        # New list of files #
        modified_files = [FASTA(self.p.modified_dir + f.filename) for f in self.refseq_bact_orig]
        # List of NCBI identifiers #
        identifiers = [g.info['Genome Name / Sample Name'][22:] for g in genomes.values() if g.fresh]
        # Filter function #
        scgc_genes_found = 0
        def only_non_scgc(reads):
            for read in reads:
                if any(i in read.description for i in identifiers):
                    scgc_genes_found =+ 1
                    continue
                else: yield read
        # Main loop #
        for orig, modif in zip(self.refseq_bact_orig, modified_files): modif.write(only_non_scgc(orig))
        # Return #
        print "Excluded %i SCGC genes" % scgc_genes_found
        self.timer.print_elapsed()
        return modified_files

    @property_cached
    def blast_db(self):
        """A blastable database of refseq + all marine organism genes"""
        blast_db = BLASTdb(self.p.genes, 'prot')
        if not self.p.genes.exists:
            # We are going to cat a whole of files together #
            print "Regrouping all fasta files together..."
            all_genes = self.refseq_bact_mod + self.refseq_arch_orig + self.marine_genomes
            shell_output("zcat %s > %s" % (' '.join(all_genes), self.p.genes))
            self.timer.print_elapsed()
            # Check that all files ended with a newline #
            print "Checking that sequence counts match..."
            assert len(blast_db) == sum(map(len,self.refseq_arch_orig)) \
                                  + sum(map(len,self.refseq_arch_orig)) \
                                  + sum(map(len,self.marine_genomes))
            self.timer.print_elapsed()
        if not self.p.pin.exists:
            # Call make DB #
            print "Building a BLAST database..."
            blast_db.makeblastdb(logfile=self.p.log)
            self.timer.print_elapsed()
        return blast_db
