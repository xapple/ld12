# Built-in modules #

# Internal modules #
from ld12.gene import Gene

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached, property_pickled
from plumbing.autopaths import AutoPaths

# Third party modules #
from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'

###############################################################################
class Genome(FASTA):
    """A FASTA file somewhere on the file system representing a genome."""

    all_paths = """
    /in_refseq_bact.pickle
    """

    def __init__(self, path, base_dir=None):
        # Attributes #
        self.path = path
        self.name = self.short_prefix
        # Extras #
        self.info = None # Filled in by the __init__.py
        self.family = None # Filled in by the family.py
        # Base directory #
        if base_dir is None: self.base_dir = self.directory + self.short_prefix + '/'
        else:                self.base_dir = base_dir
        # Auto paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def label(self):
        """This will return names such as 'HTCC1062' instead of '637000058'"""
        return self.info['taxon']

    @property_cached
    def genes(self):
        return dict((seq.id, Gene(seq, self)) for seq in self)

    @property
    def partial(self):
        """Apparently some of them are SAGs and thus only partial."""
        return True if self.filename.startswith('2236') else False

    @property
    def long_name(self):
        """A more descriptive name"""
        return "Genome '%s' (%s) of family '%s'" % (self.name, self.label, self.info['group'])

    @property
    def fresh(self):
        """Is this organism a freshwater orgiansm?"""
        return True if self.info['mf'] == 'f' else False

    @property
    def marine(self):
        """Is this organism a marine orgiansm?"""
        return not self.fresh

    @property
    def environ(self):
        """Is it fresh or marine?"""
        assert self.fresh is not self.marine
        return "fresh" if self.fresh else "marine"

    @property_pickled
    def in_refseq_bact(self):
        """Has this organism been included in the latest version of
        the refseq bacterial protein database?"""
        handle = Entrez.esearch(db="protein", retmax=10, term=self.info['taxon'])
        results = Entrez.read(handle)
        handle.close()
        count = int(results['Count'])
        return count > 0