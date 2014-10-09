# Futures #
from __future__ import division

# Built-in modules #
import os, multiprocessing
from collections import defaultdict

# Internal modules #
from ld12 import genomes
from ld12.cluster import Cluster

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.common import natural_sort, which, split_thousands
from plumbing.timer import Timer
from seqsearch.parallel import ParallelSeqSearch
from seqsearch.blast import BLASTdb

# Third party modules #
import sh, pandas
from shell_command import shell_output

###############################################################################
class Analysis(object):
    """The main object. The only mandatory argument is the input fasta file path.
    Look at the Omnigraffle outline for more information."""

    all_paths = """
    /blast/all_genes.fasta
    /blast/all_genes.fasta.nin
    /blast/all_genes.fasta.pin
    /blast/all_genes.blastout
    /mcl/bit_scores.abc
    /mcl/network.mci
    /mcl/dictionary.tab
    /mcl/clusters.txt
    /user_outputs/count_table.tsv
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
        # Paths #
        self.base_dir = out_dir
        if not self.base_dir.endswith('/'): self.base_dir += '/'
        if not os.path.exists(self.base_dir): os.makedirs(self.base_dir)
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
            shell_output('gunzip -c %s > %s' % (' '.join(genomes.values()), db))
            assert len(db) == sum(map(len,genomes.values()))
            db.makeblastdb()
            self.timer.print_elapsed()
        return db

    @property_cached
    def search(self):
        """The sequence similarity search to be run"""
        return ParallelSeqSearch(
              algorithm   = "blast",
              input_fasta = self.blast_db,
              seq_type    = self.seq_type,
              database    = self.blast_db,
              num_threads = self.num_threads,
              filtering   = {'e_value':      self.e_value,
                             'min_identity': self.min_identity,
                             'min_coverage': self.min_coverage},
              params      = {'-outfmt' : "6 qseqid sseqid bitscore pident qcovs"})

    @property
    def search_results(self):
        """For every gene, search against a database of all gene return the best hits
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

    @property_cached
    def scores(self):
        """For every gene pair that had a significant hit, what is the score of this hit."""
        result = {}
        for qseqid, sseqid, bitscore, pident, qcovs in self.search.results:
            result[(qseqid, sseqid)] = bitscore
        return result

    @property_cached
    def clusters(self):
        """A list of Clusters. See http://bioops.info/2011/03/mcl-a-cluster-algorithm-for-graphs/"""
        if not self.p.clusters.exists:
            print "Using results from %i hits" % split_thousands(len(self.scores))
            print "--> STEP 4: Running the MCL clustering"
            self.p.bit_scores.writelines(k[0]+'\t'+k[1]+'\t'+v+'\n' for k,v in self.scores.items())
            sh.mcxload("-abc", self.p.bit_scores, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            mcl = sh.Command(which('mcl'))
            mcl(self.p.network, "-I", str(self.mcl_factor), "-use-tab", self.p.dictionary, "-o", self.p.clusters)
            print "Got %i clusters" % len(self.p.clusters)
            self.timer.print_elapsed()
        # Make the clusters #
        return [Cluster(i, line, self) for i, line in enumerate(self.p.clusters)]

    @property_cached
    def best_clusters(self):
        """Subset of self.clusters. We want to find the clusters that have exactly one
        member in each of the genomes. Some genomes are partial so we will be more
        flexible on those ones."""
        self.clusters = sorted(self.clusters, key=lambda x: x.score, reverse=True)
        return self.clusters[0:50]

    #-------------------------------------------------------------------------#
    def make_trees(self):
        for c in self.best_clusters:
            print "Building tree for cluster '%s'..." % c.name
            print c.tree
            self.timer.print_elapsed()

    def save_count_table(self):
        self.count_table.to_csv(str(self.p.tsv), sep='\t', encoding='utf-8')

    @property_cached
    def count_table(self):
        """Return a dataframe with genomes as columns and clusters as rows"""
        result = defaultdict(lambda: defaultdict(int))
        for cluster in self.clusters:
            for gene in cluster.genes:
                result[cluster.name][gene.genome.name] += 1
        result = pandas.DataFrame(result, dtype=int)
        result = result.reindex_axis(sorted(result.index, key=natural_sort))
        result = result.reindex_axis(sorted(result.columns, key=natural_sort), axis='columns')
        result = result.fillna(0)
        return result

    @property
    def percent_filtered(self):
        """How many hits did we filter away ?"""
        before = sum([q.out_path.count for q in self.search.blast_queries])
        after = self.search.out_path.count
        return "%.1f%%" % (100 - (after/before)*100)
