# Built-in modules #
import re

# First party modules #

# Third party modules #

###############################################################################
class Gene(object):
    """A DNA sequence with an ID associated and belonging to a genome."""

    def __repr__(self): return '<%s object %s>' % (self.__class__.__name__, self.name)
    def __str__(self): return str(self.seq.seq)

    def __init__(self, seq, genome):
        self.seq = seq
        self.name = seq.id
        self.genome = genome
        self.annotation = None # Filled in by the __init__.py

    @property
    def long_name(self):
        """A more descriptive name"""
        return self.name + " (from " + self.genome.long_name + ")"

    @property
    def ribo_group(self):
        """If it is a ribosomal protein, what group is it part of ?"""
        results = re.findall("ribosomal protein ([LS][1-9]+)", self.annotation)
        if not results: return False
        else: return results[0]