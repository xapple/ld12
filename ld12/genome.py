# Built-in modules #

# First party modules #
from fasta import FASTA

# Third party modules #

###############################################################################
class Genome(FASTA):
    """A FASTA file somewhere on the file system."""

    @property
    def partial(self):
        """Apparently some of them are SAGs and thus only partial."""
        return True if self.filename.startswith('2236') else False