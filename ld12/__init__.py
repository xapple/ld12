b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.0.1'

# Built-in modules #
import sys, os, glob

# Internal modules #
from ld12.genome import Genome
from ld12.family import Family

# Third party modules #
import pandas

###############################################################################
# Find the data dir #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__) + '/'
repos_dir = os.path.abspath(module_dir + '../') + '/'
data_dir = module_dir + 'data/'

# Make the genomes #
fasta_paths = glob.glob(data_dir + 'genomes/*.fasta.gz')
genomes = [Genome(f) for f in fasta_paths]
genomes = dict((G.name, G) for G in genomes)

# Iterate over all genes #
genes = dict((g.name, g) for G in genomes.values() for g in G.genes.values())

# Add the metadata #
metadata = pandas.io.parsers.read_csv(data_dir + 'metadata.tsv', sep='\t', index_col=0, encoding='utf-8')
for g in genomes.values(): g.info = metadata.loc[int(g.short_prefix)]

# Add the annotations #
for line in open(data_dir + 'annotations.tsv'):
    gene_id, annotation = line.strip('\n').split('\t')
    genes[gene_id].annotation = annotation

# Make the genome families #
families = set([G.info['group'] for G in genomes.values()])
families = [Family(f, [G for G in genomes.values() if G.info['group'] == f]) for f in families]
families = dict((f.name, f) for f in families)