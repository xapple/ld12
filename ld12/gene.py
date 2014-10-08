# Built-in modules #

# First party modules #

# Third party modules #

###############################################################################
class Gene(object):
    """A sequence with an ID associated."""

    def __init__(self, seq):
        self.seq = seq
        self.id = seq.id
        self.annotation = None # Filled in by the __init__.py