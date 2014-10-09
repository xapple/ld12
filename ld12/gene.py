# Built-in modules #

# First party modules #

# Third party modules #

###############################################################################
class Gene(object):
    """A sequence with an ID associated."""

    def __repr__(self): return '<%s object %s>' % (self.__class__.__name__, self.id)
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
