# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #
from ld12.cluster import Cluster
from ld12.genome import Genome

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.common import natural_sort, which
from parallelblast import BLASTdb, BLASTquery
from parallelblast.results import tabular_keys

# Third party modules #
import sh, pandas
from shell_command import shell_output
from tqdm import tqdm

###############################################################################
class Analysis(object):
    """All sequences from all genomes are blasted against themselves.
    Then we use the MCL algorithm (Markov Cluster Algorithm) to form clusters.
    Then we count the genomes with the right number of single copy genes."""

    blast_params = {'-e': 0.001, '-m': 8}
    minimum_identity = 30.0
    mimimum_coverage = 50.0
    sequence_type ='aminoacid' or 'nucleotide'

    all_paths = """
    /all_sequences.fasta
    /all_sequences.fasta.nin
    /all_sequences.fasta.pin
    /all_sequences.blastout
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

    def __init__(self, genomes, base_dir='.'):
        # Attributes #
        self.genomes = genomes
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def blast_db(self):
        """A blastable database of all genes"""
        assert self.genomes
        dbtype = 'nucl' if self.sequence_type == 'nucleotide' else 'prot'
        db = BLASTdb(self.p.all_fasta, dbtype=dbtype)
        if not self.p.all_nin and not self.p.all_pin:
            print "Building BLASTable database with all genes..."
            shell_output('cat %s > %s' % (' '.join(self.genomes), db))
            assert len(db.ids) == len(set(db.ids))
            db.makeblastdb()
        return db

    @property_cached
    def query(self):
        """The blast query to be executed"""
        algorithm = 'blastn' if self.sequence_type == 'nucleotide' else 'blastp'
        return BLASTquery(self.blast_db, self.blast_db, self.blast_params, algorithm, 'legacy')

    @property_cached
    def blastout(self):
        """The blast results"""
        if not self.p.all_blastout:
            print "Self-BLASTing database '%s'..." % self.blast_db.relative_path
            self.query.run()
        return self.p.all_blastout

    @property_cached
    def filtered(self):
        """We want to check the percent identify and the coverage of the hit to the query"""
        def good_iterator(blastout):
            print "Filtering BLAST hits..."
            for line in tqdm(blastout, total=len(blastout)):
                info = dict(zip(tabular_keys, line.split()))
                if float(info['perc_identity']) < self.minimum_identity: continue
                query_cov = (float(info['query_end']) - float(info['query_start']))
                query_cov = 100.0 * abs(query_cov / self.blast_db.length_by_id[info['query_id']])
                subj_cov  = (float(info['subject_end']) - float(info['subject_start']))
                subj_cov  = 100.0 * abs(subj_cov  / self.blast_db.length_by_id[info['subject_id']])
                coverage = min(query_cov, subj_cov)
                if coverage < self.mimimum_coverage: continue
                yield line
        if not self.p.filtered_blastout:
            print "Making SQLite database with reads from '%s'..." % self.blastout.relative_path
            print "Result in '%s'." % self.blast_db.sql.relative_path
            self.p.filtered_blastout.writelines(good_iterator(self.blastout))
        return self.p.filtered_blastout

    @property
    def percent_filtered(self):
        """How many hits did we filter away ?"""
        percentage = lambda x,y: (len(x)/len(y))*100 if len(y) != 0 else 0
        return "%.1f%%" % (100 - percentage(self.filtered, self.blastout))

    @property_cached
    def clusters(self):
        """A list of Clusters. See http://bioops.info/2011/03/mcl-a-cluster-algorithm-for-graphs/"""
        if not self.p.clusters.exists:
            print "Running the MCL clustering on '%s'..." % self.filtered.relative_path
            shell_output("cut -f 1,2,11 %s > %s" % (self.filtered, self.p.filtered_abc))
            sh.mcxload("-abc", self.p.filtered_abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            mcl = sh.Command(which('mcl'))
            mcl(self.p.network, "-I", "1.5", "-use-tab", self.p.dictionary, "-o", self.p.clusters)
        return [Cluster(i, frozenset(line.split()), self) for i, line in enumerate(self.p.clusters)]

    @property_cached
    def count_table(self):
        """Genomes as columns and clusters as rows"""
        result = defaultdict(lambda: defaultdict(int))
        for genome in self.genomes:
            for cluster in self.clusters:
                for id_num in cluster.ids:
                    if id_num in genome:
                        result[genome.prefix][cluster.name] += 1
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