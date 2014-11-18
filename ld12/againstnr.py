# Built-in modules #

# Internal modules #
import socket, os

# First party modules #
from seqsearch.parallel import ParallelSeqSearch
from seqsearch.blast import BLASTdb
from plumbing.cache import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'
host = socket.gethostname()

# Hardcoded database location #
if host.startswith('milou'): refseq_special_db = "/gulo/glob/alexe/databases/refseq/"
else:                        refseq_special_db = home + "LD12/databases/refseq"

###############################################################################
class AgainstNR(object):
    """This sub-object takes care of BLASTing the genes against the NR database
    and parsing the results form that. We want to blast only the fresh water
    clusters against a modified refseq database."""

    all_paths = """
    /blast/fresh_genes.fasta
    /blast/fresh_genes.blastout
    """

    def __init__(self, analysis,
                 e_value      = 0.001,
                 min_identity = 0.3,
                 min_coverage = 0.5):
        # Attributes #
        self.analysis = analysis
        self.clusters = analysis.fresh_clusters
        self.seq_type = analysis.seq_type
        # Other #
        self.e_value      = e_value
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        # Paths #
        self.base_dir = analysis.p.againstnr_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def fresh_fasta(self):
        """A file containing all the fresh water genes"""
        return BLASTdb(refseq_special_db, self.seq_type)

    @property_cached
    def blast_db(self):
        """A blastable database of all refseq + all marine organism genes"""
        return BLASTdb(refseq_special_db, self.seq_type)

    @property_cached
    def search(self):
        """The sequence similarity search to be run"""
        return ParallelSeqSearch(
              algorithm   = "blast",
              input_fasta = self.fresh_fasta,
              seq_type    = self.seq_type,
              database    = self.blast_db,
              num_threads = self.num_threads,
              filtering   = {'e_value':      self.e_value,
                             'min_identity': self.min_identity,
                             'min_coverage': self.min_coverage},
              params      = {'-outfmt' : "6 qseqid sseqid bitscore pident qcovs"})

    @property
    def search_results(self):
        """For every gene, search against a database of all gene, return the best hits
        after filtering."""
        # Check that the search was run #
        if not self.search.out_path.exists:
            print "Using: %i genes" % split_thousands(len(self.blast_db))
            print "--> STEP 2: Similarity search against all genes with %i processes" % self.num_threads
            self.search.run()
            self.timer.print_elapsed()
            print "--> STEP 3: Filter out bad hits from the search results"
            self.search.filter()
            if self.search.out_path.count_bytes == 0:
                raise Exception("Found exactly zero hits after the similarity search.")
            print "Filtered %s of the hits" % self.percent_filtered
            self.timer.print_elapsed()
        # Parse the results #
        return self.search.results
