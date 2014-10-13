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
        # Basic #
        self.name = name
        self.genomes = genomes
        # Add the link from genomes to Family #
        for genome in self.genomes: genome.family = self

    @property
    def genes(self):
        """All the genes of all the genomes in this family"""
        return [g for genome in self.genomes for g in genome.genes.values()]
