# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.cache import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class AgainstNR(object):
    """This sub-object takes care of BLASTing the genes against the NR database
    and parsing the results form that."""

    def __init__(self, analysis):
        # Links
        self.analysis = analysis
