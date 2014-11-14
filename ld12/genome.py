# Built-in modules #

# Internal modules #
from ld12.gene import Gene

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class Genome(FASTA):
    """A FASTA file somewhere on the file system representing a genome."""

    def __init__(self, path):
        self.path = path
        self.name = self.short_prefix
        self.info = None # Filled in by the __init__.py
        self.family = None # Filled in by the family.py

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
        return "Genome '" + self.label + "'' of family '" + self.info['group'] + "'"

    @property
    def fresh(self):
        """Is this organism a freshwater orgiansm."""
        return True if self.info['mf'] == 'f' else False
