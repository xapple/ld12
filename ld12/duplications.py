# Built-in modules #
import socket, os

# Internal modules #
from ld12 import genomes, genes

# First party modules #
from seqsearch.parallel import ParallelSeqSearch
from seqsearch.blast import BLASTdb
from seqsearch.common import UtilsNCBI
from plumbing.cache import property_cached
from plumbing.autopaths import AutoPaths
from plumbing.common import split_thousands
from fasta import FASTA

# Third party modules #
from shell_command import shell_output
from tqdm import tqdm

# Constants #
home = os.environ['HOME'] + '/'
host = socket.gethostname()

# Hardcoded database location #
if host.startswith('m'): refseq_special_db = "/gulo/glob/alexe/databases/refseq/refseq_SAR11"
else:                    refseq_special_db = home + "LD12/databases/refseq/refseq_SAR11"

###############################################################################
class Duplications(object):
    """This sub-object takes care of BLASTing the genes against the refseq database
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
        self.analysis    = analysis
        self.clusters    = analysis.fresh_clusters
        self.seq_type    = analysis.seq_type
        self.num_threads = analysis.num_threads
        self.timer       = analysis.timer
        # Other #
        self.e_value      = e_value
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        # Paths #
        self.base_dir = analysis.p.duplications_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def fresh_fasta(self):
        """A file containing all the fresh water genes"""
        fasta = FASTA(self.p.fresh_fasta)
        if not fasta.exists:
            print "Building fasta file with all fresh genes..."
            fresh = [g for g in genomes.values() if g.fresh]
            shell_output('gunzip -c %s > %s' % (' '.join(fresh), fasta))
            assert len(fasta) == sum(map(len, fresh))
            self.timer.print_elapsed()
        return fasta

    @property_cached
    def blast_db(self):
        """A blastable database of all refseq + all marine organism genes"""
        return BLASTdb(refseq_special_db, 'nucl')

    @property_cached
    def search(self):
        """The sequence similarity search to be run (tblastn)"""
        return ParallelSeqSearch(
              algorithm   = "blast",
              input_fasta = self.fresh_fasta,
              seq_type    = self.seq_type,
              database    = self.blast_db,
              num_threads = self.num_threads,
              filtering   = {'max_targets':  1,
                             'e_value':      self.e_value,
                             'min_identity': self.min_identity,
                             'min_coverage': self.min_coverage},
              params      = {'-outfmt' : "6 qseqid sseqid bitscore pident qcovs"})

    @property
    def search_results(self):
        """For every gene, search against a database of all gene, return the best hit
        after filtering."""
        # Check that the search was run #
        if not self.search.out_path.exists:
            print "Using: %s genes" % split_thousands(len(self.fresh_fasta))
            print "Similarity search against custom database for all fresh genes with %i processes" % self.num_threads
            self.search.run()
            self.timer.print_elapsed()
            print "Filter out bad hits from the search results"
            self.search.filter()
            if self.search.out_path.count_bytes == 0:
                raise Exception("Found exactly zero hits after the similarity search.")
            self.timer.print_elapsed()
        # Parse the results #
        return self.search.results

    def assign_best_hits(self):
        """Parse the results and add the best hit information for each Gene
        object in each freshwater Genome object"""
        last_query_id = -1
        for query_id, hit_id, bitscore, identity, coverage in self.search_results:
            if query_id != last_query_id:
                gene = genes[query_id]
                gene.best_hit = hit_id
                last_query_id = query_id
                continue

    def assign_taxonomy(self):
        """Use the best hit information for each Gene object to add the assign_taxonomy
        information of each best hit to each Gene object"""
        # Check all genes have a hit #
        self.fresh_genes      = [g for g in genes.values() if g.genome.fresh]
        self.hit_genes        = [g for g in self.fresh_genes if hasattr(g, 'best_hit')]
        self.no_hit_genes     = [g for g in self.fresh_genes if not hasattr(g, 'best_hit')]
        self.ncbi_hit_genes   = [g for g in self.hit_genes if g.best_hit.startswith('gi')]
        self.marine_hit_genes = [g for g in self.hit_genes if g.best_hit.startswith('Pelub58DRAFT')]
        # Check there are no others #
        assert sum(map(len, (self.marine_hit_genes, self.ncbi_hit_genes, self.no_hit_genes))) == len(self.fresh_genes)
        # Print those that don't #
        for gene in self.no_hit_genes: print "Gene %s did not get a best hit against Refseq+Marine" % gene.name
        # Assign taxonomy #
        for gene in self.marine_hit_genes: gene.taxonomy = "Pelub58"
        # Download in batch #
        for gene in self.ncbi_hit_genes: gene.gi_num = gene.best_hit.split('|')[1]
        ncbi = UtilsNCBI()
        taxonomies = ncbi.gi_num_to_tax([g.gi_num for g in self.ncbi_hit_genes])
        for gene in self.ncbi_hit_genes: gene.taxonomy = taxonomies[gene.gi_num]

    def save_taxonomy_info(self):
        """Save the results"""
        pass

    def make_plot(self):
        """Plot the hit information"""
        pass