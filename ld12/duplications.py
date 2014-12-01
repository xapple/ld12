# Built-in modules #
import os
from collections import OrderedDict

# Internal modules #
from ld12 import genomes, genes, families
from ld12.refseq import RefSeqProkPlusMarine

# First party modules #
from seqsearch.parallel import ParallelSeqSearch
from seqsearch.common import UtilsNCBI
from plumbing.cache import property_cached, property_pickled
from plumbing.autopaths import AutoPaths
from plumbing.common import split_thousands
from plumbing.graphs import Graph
from fasta import FASTA

# Third party modules #
import pandas
from shell_command import shell_output

# Plot #
import matplotlib
matplotlib.use('Agg', warn=False)
from matplotlib import pyplot

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Duplications(object):
    """This sub-object takes care of BLASTing the genes against the refseq database
    and parsing the results form that. We want to blast only the fresh water
    clusters against a modified refseq database."""

    all_paths = """
    /blast/fresh_genes.fasta
    /blast/fresh_genes.blastout
    /ncbi_taxonomy/gi_to_record.pickle
    /database/
    /graphs/
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
        # The database #
        self.refseq = RefSeqProkPlusMarine(self.p.database_dir, self)
        # The final plot #
        self.plot = TaxonomyPlot(self)

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
    def search(self):
        """The sequence similarity search to be run"""
        return ParallelSeqSearch(
            algorithm    = "blast",
            input_fasta  = self.fresh_fasta,
            seq_type     = self.seq_type,
            database     = self.refseq.blast_db,
            num_parts    = 200,
            filtering    = {'max_targets' : 50,
                            'e_value'     : self.e_value,
                            'min_identity': self.min_identity,
                            'min_coverage': self.min_coverage},
            params       = {'-outfmt'     : '"6 qseqid sseqid bitscore pident qcovs"'},
            slurm_params = {'time'        : '7-00:00:00',
                            'cores'       : 1,
                            'project'     : 'b2011035',
                            'partition'   : 'core'})

    @property
    def search_results(self):
        """For every gene, search against a database of all gene, return the best hit
        after filtering."""
        # Check that the search was run #
        if not self.search.out_path.exists:
            print "Using: %s genes" % split_thousands(len(self.fresh_fasta))
            print "Similarity search against custom database for all fresh genes with %i processes" % self.num_threads
            self.search.run_local()
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
        # Message #
        if self.search_results: print "Parsing the best hits file..."
        # Only one best hit per gene #
        last_query_id = -1
        for query_id, hit_id, bitscore, identity, coverage in self.search_results:
            if query_id != last_query_id:
                gene = genes[query_id]
                gene.best_hit = hit_id
                last_query_id = query_id
                continue
        # Make groups #
        self.fresh_genes      = [g for g in genes.values() if g.genome.fresh]
        self.hit_genes        = [g for g in self.fresh_genes if hasattr(g, 'best_hit')]
        self.no_hit_genes     = [g for g in self.fresh_genes if not hasattr(g, 'best_hit')]
        self.ncbi_hit_genes   = [g for g in self.hit_genes if g.best_hit.startswith('gi')]
        self.marine_hit_genes = [g for g in self.hit_genes if not g.best_hit.startswith('gi')]
        # Check there are no others #
        assert sum(map(len, (self.marine_hit_genes, self.ncbi_hit_genes, self.no_hit_genes))) == len(self.fresh_genes)
        # Extract numbers #
        for gene in self.ncbi_hit_genes: gene.gi_num = gene.best_hit.split('|')[1]
        # Extract gene objects #
        for gene in self.marine_hit_genes: gene.marine_hit = genes[gene.best_hit]
        # All possible GI numbers #
        self.all_gi_nums = [g.gi_num for g in self.ncbi_hit_genes]
        # Done #
        self.timer.print_elapsed()

    @property_pickled
    def gi_to_record(self):
        """Link all possible GIs we found to their record and save the result"""
        # Download in batch #
        ncbi = UtilsNCBI()
        return ncbi.gis_to_records([g.gi_num for g in self.ncbi_hit_genes])

    def assign_taxonomy(self):
        """Use the best hit information for each Gene object to add the taxonomy
        information of each best hit to each Gene object"""
        print "Linking GI numbers to their NCBI taxonomy..."
        ncbi = UtilsNCBI()
        for g in self.no_hit_genes:     g.taxonomy = None
        for g in self.marine_hit_genes: g.taxonomy = g.marine_hit.genome.family.name
        for g in self.ncbi_hit_genes:   g.taxonomy = ncbi.record_to_taxonomy(self.gi_to_record[g.gi_num])
        self.timer.print_elapsed()

    @property
    def duplications_stats(self):
        yield """
        Total fresh water genes: %s
        Did not get a top hit: %s
        Did get a top hit against one of the marine genes: %s
        Did get a top hit against one of the other NCBI genes: %s
        \n\n""" % (len(self.fresh_genes),
                   len(self.no_hit_genes),
                   len(self.marine_hit_genes),
                   len(self.ncbi_hit_genes))
        yield "Gene name\tHit type\tHit reference\tHit taxonomy\n"
        for g in self.no_hit_genes:     yield g.name + "\t" + "No hit" + '\t' + "None" + "\t" + "None" + '\n'
        for g in self.marine_hit_genes: yield "\t".join((g.name, "Marine hit", g.marine_hit.name, g.taxonomy)) + '\n'
        for g in self.ncbi_hit_genes:   yield "\t".join((g.name, "NCBI hit",   g.gi_num,          g.taxonomy)) + '\n'
        yield '\n'

    def save_duplications_stats(self):
        """Save the results"""
        self.analysis.p.duplications.writelines(self.duplications_stats)

###############################################################################
class TaxonomyPlot(Graph):
    short_name = 'fresh_taxonomy'
    bottom = 0.35
    left   = 0.1

    def plot(self):
        # First the no hits #
        no_hits = {"No hits": len(self.parent.no_hit_genes)}
        # The marine hits #
        fams = OrderedDict([(f,0) for f in families])
        for g in self.parent.marine_hit_genes: fams[g.taxonomy] += 1
        # Then the ncbi hits #
        categories = OrderedDict((('Life',0),
                                  ('Bacteria',0),
                                  ('Proteobacteria',0),
                                  ('Alphaproteobacteria',0),
                                  ('SAR11 cluster',0),
                                  ('Candidatus Pelagibacter',0)))
        for g in self.parent.ncbi_hit_genes:
            tax = g.taxonomy.split(';')
            tax.append("")
            tax.append("")
            tax.append("")
            if   tax[0].strip() != 'Bacteria':                 categories['Life']                     += 1
            elif tax[1].strip() != 'Proteobacteria':           categories['Bacteria']                 += 1
            elif tax[2].strip() != 'Alphaproteobacteria':      categories['Proteobacteria']           += 1
            elif tax[3].strip() != 'SAR11 cluster':            categories['Alphaproteobacteria']      += 1
            elif tax[4].strip() != 'Candidatus Pelagibacter':  categories['SAR11 cluster']            += 1
            else:                                              categories['Candidatus Pelagibacter']  += 1
        # Frame #
        self.frame = OrderedDict()
        self.frame.update(no_hits)
        self.frame.update(categories)
        self.frame.update(fams)
        self.frame = pandas.Series(self.frame)
        # Plot #
        axes = self.frame.plot(kind='bar', color='gray')
        fig = pyplot.gcf()
        axes.set_title("Taxonomy distribution for the best hit against refseq for all freshwater genes")
        axes.set_ylabel("Number of best hits with this taxonomy")
        axes.xaxis.grid(True)
        axes.set_yscale('symlog')
        self.save_plot(fig, axes, sep=('y',))
        pyplot.close(fig)