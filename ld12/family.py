# Built-in modules #

# Internal modules #

# First party modules #

# Third party modules #

###############################################################################
class Family(object):
    """A collection of genomes that all share the same family."""

    def __repr__(self): return '<%s object "%s" with %i genomes>' % \
        (self.__class__.__name__, self.name, len(self.genomes))

    def __init__(self, name, genomes):
        self.name = name
        self.genomes = genomes