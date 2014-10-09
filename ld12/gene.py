# Built-in modules #

# First party modules #

# Third party modules #

###############################################################################
class Gene(object):
    """A sequence with an ID associated."""

    def __repr__(self): return '<%s object %s>' % (self.__class__.__name__, self.id)

    def __init__(self, seq, genome):
        self.seq = seq
        self.name = seq.id
        self.genome = genome
        self.annotation = None # Filled in by the __init__.py