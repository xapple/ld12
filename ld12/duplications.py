# Built-in modules #
import os, re
from collections import OrderedDict

# Plot #
import matplotlib
matplotlib.use('Agg', warn=False)
from matplotlib import pyplot

# Internal modules #
from ld12 import genomes, genes, families
from ld12.refseq import RefSeqProkPlusMarine

# First party modules #
from seqsearch.parallel import ParallelSeqSearch
from seqsearch.ncbi import UtilsNCBI
from plumbing.cache import property_cached, property_pickled
from plumbing.autopaths import AutoPaths
from plumbing.common import split_thousands
from plumbing.graphs import Graph
from fasta import FASTA

# Third party modules #
import pandas
from shell_command import shell_output
from tqdm import tqdm

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Duplications(object):
    """This sub-object takes care of BLASTing the genes against the refseq database
    and parsing the results from that. We want to blast only the fresh water
    clusters against a modified refseq database."""

    all_paths = """
    /blast/fresh_genes.fasta
    /blast/fresh_genes.blastout
    /ncbi_taxonomy/gi_to_record.pickle
    /database/
    /graphs/
    /user_outputs/hit_stats.txt
    /user_outputs/duplications_stats.tsv
    /user_outputs/blast_stats.tsv
    """

    def __init__(self, analysis, genes_to_consider, base_dir,
                 e_value           = 0.001,
                 min_identity      = 0.3,
                 min_coverage      = 0.5):
        # Attributes #
        self.analysis    = analysis
        self.seq_type    = analysis.seq_type
        self.num_threads = analysis.num_threads
        self.timer       = analysis.timer
        # Other #
        self.e_value      = e_value
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        # Paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # The database #
        self.refseq = RefSeqProkPlusMarine(self.p.database_dir, self)
        # The genes #
        self.genes = genes_to_consider
        self.genome_names = [G.info['taxon'] for G in genomes.values()]
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
        """The sequence similarity search to be run."""
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
        """Return the best hits after filtering."""
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

    @property_cached
    def all_gis(self):
        """Extract all GI numbers we found in the search results"""
        # Message #
        if self.search_results: print "Parsing the best hits file..."
        # Main loop #
        gis = set()
        for line in self.search_results:
            hit_id = line[1]
            if not hit_id.startswith('gi'): continue
            gi_num = hit_id.split('|')[1]
            gis.add(gi_num)
        return gis
        # Done #
        self.timer.print_elapsed()

    @property_pickled
    def gi_to_record(self):
        """Link all possible GIs we found to their records and save the result."""
        # Message #
        print "Downloading records from NCBI for %i GIs..." % len(self.all_gis)
        # Download in batch #
        ncbi = UtilsNCBI()
        gis = list(self.all_gis)
        return ncbi.gis_to_records(gis)
        # Done #
        self.timer.print_elapsed()

    #-------------------------------------------------------------------------#
    #                                 HITS                                    #
    #-------------------------------------------------------------------------#
    def assign_hits(self):
        """Parse the results and add the hits information for each Gene
        object in each freshwater Genome object."""
        # Just assign #
        for query_id, hit_id, bitscore, identity, coverage in self.search_results:
            gene = genes[query_id].raw_hits.append((hit_id, bitscore, coverage))
        # Load the pickled property #
        print "Loading all NCBI records into RAM..."
        print "Got %i records." % len(self.gi_to_record)
        self.timer.print_elapsed()
        # Process #
        print "Processing top hits for %i genes..." % len(self.genes)
        for gene in tqdm(self.genes):
            result = []
            for hit_id, bitscore, coverage in gene.raw_hits:
                # The id and score #
                hit = {'id':       hit_id,
                       'score':    bitscore,
                       'coverage': coverage}
                # The source: refseq, missing #
                if hit_id.startswith('gi'):
                    hit['source'] = "refseq"
                else:
                    hit['source'] = "missing"
                    hit['gene']   = genes[hit_id]
                    hit['genome'] = hit['gene'].genome
                # The type: marine, fresh, other #
                if hit['source'] == "missing":
                    hit['type'] = hit['gene'].genome.environ
                if hit['source'] == "refseq":
                    hit['gi_num'] = hit_id.split('|')[1]
                    hit['record'] = self.gi_to_record[hit['gi_num']]
                    hit['header'] = hit['record']['GBSeq_organism']
                    if not any(name in hit['header'] for name in self.genome_names):
                        hit['type'] = 'other'
                    else:
                        matching_genomes = [G for G in genomes.values() if re.search(G.info['taxon'] + '$', hit['header'])]
                        assert len(matching_genomes) == 1
                        hit['genome'] = matching_genomes[0]
                        hit['type']   = hit['genome'].environ
                # The taxonomy #
                if hit['type'] == "other": hit['taxonomy'] = hit['record']['GBSeq_taxonomy']
                else:                      hit['taxonomy'] = hit['genome'].family.name
                # Skip the ones that are classified as MULTISPECIES #
                if hit['source'] == "refseq" and 'MULTISPECIES' in hit['record']['GBSeq_definition']: continue
                # Append it #
                result.append(hit)
                # Stop at the first non-fresh #
                if hit['type'] != "fresh": break
            gene.hits = result
        # Group into categories #
        self.no_hit_genes   = [g for g in self.genes if len(g.raw_hits) == 0]
        self.yes_hit_genes  = [g for g in self.genes if len(g.raw_hits) != 0]
        self.top_is_fresh   = [g for g in self.yes_hit_genes if g.hits[0]['type'] == 'fresh']
        self.best_is_marine = [g for g in self.yes_hit_genes if g.hits[-1]['type'] == 'marine']
        self.best_is_other  = [g for g in self.yes_hit_genes if g.hits[-1]['type'] == 'other']
        # Assign the first deviating taxa #
        for g in self.best_is_other:
            tax = g.hits[-1]['taxonomy'].split(';')
            tax += ["", "", ""]
            if   tax[0].strip() != 'Bacteria':                 g.best_tax = 'Life'
            elif tax[1].strip() != 'Proteobacteria':           g.best_tax = 'Bacteria'
            elif tax[2].strip() != 'Alphaproteobacteria':      g.best_tax = 'Proteobacteria'
            elif tax[3].strip() != 'SAR11 cluster':            g.best_tax = 'Alphaproteobacteria'
            elif tax[4].strip() != 'Candidatus Pelagibacter':  g.best_tax = 'SAR11 cluster'
            else:                                              g.best_tax = 'Candidatus Pelagibacter'
        # Done #
        self.timer.print_elapsed()

    @property
    def hit_stats(self):
        yield """
        Total fresh water genes: %s
        Did not get any hits all: %s
        The absolute top hit is against a fresh water genome: %s
        The first non-fresh hit is against one of the marine genomes: %s
        The first non-fresh hit is against one of the other NCBI genomes: %s
        \n\n""" % (len(self.genes),
                   len(self.no_hit_genes),
                   len(self.top_is_fresh),
                   len(self.best_is_marine),
                   len(self.best_is_other))

    def save_hit_stats(self):
        """Save the results"""
        self.p.hit_stats.writelines(self.hit_stats)

    #-------------------------------------------------------------------------#
    #                             DUPLICATIONS                                #
    #-------------------------------------------------------------------------#
    @property
    def duplications_stats(self):
        result = OrderedDict()
        for g in self.genes:
            # Basic stats #
            result[g.name] = OrderedDict()
            result[g.name]['genome']                = g.genome.name
            result[g.name]['taxon']                 = g.genome.info['taxon']
            result[g.name]['# of hits']             = len(g.raw_hits)
            result[g.name]['# of fresh hits']       = len([h for h in g.hits if h['type'] == 'fresh'])
            result[g.name]['Is there a marine hit'] = len([h for h in g.hits if h['type'] == 'marine'])
            result[g.name]['Is there a refseq hit'] = len([h for h in g.hits if h['type'] == 'other'])
            # The hits not in this genome #
            fresh_outsiders = [h for h in g.hits if h['type'] == 'fresh' and h['genome'] is not g.genome]
            result[g.name]['# of fresh hits not in genome'] = len(fresh_outsiders)
        # Make a dataframe #
        result = pandas.DataFrame.from_dict(result)
        result = result.transpose()
        return result

    def save_duplications_stats(self):
        """Save the dataframe above in a CSV file"""
        self.duplications_stats.to_csv(str(self.p.duplications), sep='\t', encoding='utf-8')

    #-------------------------------------------------------------------------#
    @property
    def blast_stats(self):
        columns = ['Query gene', 'Query genome', 'Query taxon',
                   'Query length', 'Query coverage',
                   'Hit number', 'Hit Type', 'Hit source', 'Hit ID', 'Hit bit score',
                   'Hit Taxonomy', 'Hit is in same genome as query']
        result = []
        # Main loop #
        for g in self.genes:
            if len(g.hits) == 0:
                result.append((g.name, g.genome.name, g.genome.info['taxon'],
                               len(g.seq), "nohit",
                               -1 , "nohit", "nohit", "nohit",  "nohit",
                               "nohit", "nohit"))
                continue
            for i, hit in enumerate(g.hits):
                result.append((g.name, g.genome.name, g.genome.info['taxon'],
                               len(g.seq), hit['coverage'],
                               i+1 , hit['type'], hit['source'], hit['id'], hit['score'],
                               hit['taxonomy'], hit['type'] == 'fresh' and hit['genome'] is g.genome))
        # Make a data frame #
        result = pandas.DataFrame(result, columns=columns)
        return result
        pass

    def save_blast_stats(self):
        """Save the dataframe above in a CSV file"""
        self.blast_stats.to_csv(str(self.p.blast_stats), sep='\t', encoding='utf-8')

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
        for g in self.parent.best_is_marine: fams[g.hits[-1]['taxonomy']] += 1
        # Then the ncbi hits #
        self.categories = OrderedDict((('Life',                    0),
                                       ('Bacteria',                0),
                                       ('Proteobacteria',          0),
                                       ('Alphaproteobacteria',     0),
                                       ('SAR11 cluster',           0),
                                       ('Candidatus Pelagibacter', 0)))
        for g in self.parent.best_is_other: self.categories[g.best_tax] += 1
        # Frame #
        self.frame = OrderedDict()
        self.frame.update(no_hits)
        self.frame.update(self.categories)
        self.frame.update(fams)
        self.frame = pandas.Series(self.frame)
        # Plot #
        axes = self.frame.plot(kind='bar', color='gray')
        fig = pyplot.gcf()
        axes.set_title("Taxonomy distribution for the best hit against refseq for all freshwater genes")
        axes.set_ylabel("Number of best hits with this taxonomy")
        axes.xaxis.grid(True)
        self.save_plot(fig, axes, sep=('y',))
        pyplot.close(fig)