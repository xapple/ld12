b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.0.1'

# Built-in modules #
import sys, os, glob

# Internal modules #
from ld12.genome import Genome

# Third party modules #
import pandas

# Find the data dir #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__) + '/'
repos_dir = os.path.abspath(module_dir + '../') + '/'
data_dir = module_dir + 'data/'

# Make the genomes #
fasta_paths = glob.glob(data_dir + 'genomes/*.fasta.gz')
genomes = [Genome(f) for f in fasta_paths]

# Add the metadata #
metadata = pandas.io.parsers.read_csv(data_dir + 'metadata.tsv', sep='\t', index_col=0, encoding='utf-8')
for g in genomes: g.info = metadata.loc[int(g.short_prefix)]