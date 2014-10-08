# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict
import multiprocessing

# Internal modules #
from ld12.cluster import Cluster

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.common import natural_sort, which
from plumbing.timer import Timer
from seqsearch import SeqSearch
from seqsearch.blast import BLASTdb

# Third party modules #
import sh, pandas
from shell_command import shell_output

###############################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.
    Look at the Omnigraffle outline for more information."""

    all_paths = """
    /all_genes.fasta
    /all_genes.fasta.nin
    /all_genes.fasta.pin
    /all_genes.blastout
    /filtered.blastout
    /filtered.abc
    /network.mci
    /dictionary.tab
    /clusters.txt
    /count_table.tsv
    /master.aln
    /master.tree
    /clusters/
    """

    def __repr__(self): return '<%s object with %i genomes>' % \
        (self.__class__.__name__, len(self.genomes))

    def __init__(self,
                 out_dir      = './output/',
                 e_value      = 0.001,
                 min_identity = 0.3,
                 min_coverage = 0.5,
                 seq_type     = 'prot' or 'nucl',
                 mcl_factor   = 1.5,
                 num_threads  = None):
        # Paths #
        self.base_dir = out_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Other #
        self.e_value      = e_value
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.seq_type     = seq_type
        self.mcl_factor   = mcl_factor
        # Number of cores to use #
        if num_threads is None: self.num_threads = min(multiprocessing.cpu_count(), 32)
        else: self.num_threads = int(num_threads)
        # Time the pipeline execution #
        self.timer = Timer()

    @property_cached
    def blast_db(self):
        """A blastable database of all genes"""
        db = BLASTdb(self.p.all_fasta, self.seq_type)
        if not self.p.all_nin and not self.p.all_pin:
            print "--> STEP 1: Building BLASTable database with all genes..."
            shell_output('cat %s > %s' % (' '.join(self.genomes), db))
            assert len(db.ids) == len(set(db.ids))
            db.makeblastdb()
            self.timer.print_elapsed()
        return db

    @property_cached
    def filtering(self):
        """Return a dictionary with the filtering options for the sequence similarity
        search."""
        return {'e_value':      self.e_value,
                'min_identity': self.min_identity,
                'min_coverage': self.min_coverage}

    @property_cached
    def search(self):
        """The sequence similarity search to be run"""
        return SeqSearch(input_fasta = self.blast_db,
                         seq_type    = self.seq_type,
                         database    = self.search_db,
                         algorithm   = "blast",
                         num_threads = self.num_threads,
                         filtering   = self.filtering,
                         params      = {'-outfmt' : "6 qseqid sseqid bitscore pident qcovs"})

    @property
    def search_results(self):
        """For every gene, search against a database of all gene return the best hits
        after filtering."""
        # Check that the search was run #
        if not self.search.out_path.exists:
            print "Using: %i genes" + len(self.blast_db)
            print "--> STEP 2: Similarity search against all genes"
            self.search.run()
            self.timer.print_elapsed()
            print "--> STEP 3: Filter out bad hits from the search results"
            self.search.filter()
            self.timer.print_elapsed()
            if self.search.out_path.count_bytes == 0:
                raise Exception("Found exactly zero hits after the similarity search.")
        # Parse the results #
        return self.search.results

    @property_cached
    def scores(self):
        """For every gene pair that had a significant hit, what is the score of this hit."""
        result = {}
        for line in self.search.results:
            qseqid, sseqid, bitscore, pident, qcovs = line.split()
            result[(qseqid, sseqid)] = bitscore
        return result

    @property_cached
    def clusters(self):
        """A list of Clusters. See http://bioops.info/2011/03/mcl-a-cluster-algorithm-for-graphs/"""
        if not self.p.clusters.exists:
            print "Using results from %i hits" % len(self.scores)
            print "--> STEP 4: Running the MCL clustering"
            self.p.filtered_abc.writelines(k[0]+'\t'+k[1]+'\t'+v+'\n' for k,v in self.scores)
            sh.mcxload("-abc", self.p.filtered_abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            mcl = sh.Command(which('mcl'))
            mcl(self.p.network, "-I", str(self.mcl_factor), "-use-tab", self.p.dictionary, "-o", self.p.clusters)
            self.timer.print_elapsed()
        # Make the clusters #
        return [Cluster(i, line) for i, line in enumerate(self.p.clusters)]

    @property_cached
    def count_table(self):
        """Return a dataframe with genomes as columns and clusters as rows"""
        result = defaultdict(lambda: defaultdict(int))
        for cluster in self.clusters:
            for gene in cluster.genes:
                result[cluster.name][gene.genome.prefix] += 1
        result = pandas.DataFrame(result)
        result = result.reindex_axis(sorted(result.index, key=natural_sort))
        result = result.fillna(0)
        return result

    def save_count_table(self):
        self.count_table = self.count_table.reindex([c.name for c in self.clusters])
        self.count_table.to_csv(str(self.p.tsv), sep='\t', encoding='utf-8')

    @property_cached
    def single_copy_clusters(self):
        """Subset of self.clusters. Which clusters appear exactly once in each genome.
        Some genomes are partial so we will be more flexible on those ones."""
        self.clusters = sorted(self.clusters, key=lambda x: x.score, reverse=True)
        return self.clusters[0:100]