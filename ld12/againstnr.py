# Built-in modules #

# Internal modules #
import socket, os

# First party modules #
from plumbing.cache import property_cached
#from plumbing.autopaths import AutoPaths

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'
host = socket.gethostname()

# Hardcoded database location #
if host.startswith('milou'): refseq_special_db = "/gulo/glob/alexe/databases/refseq/"
else:                        refseq_special_db = home + "LD12/databases/refseq"

###############################################################################
class AgainstNR(object):
    """This sub-object takes care of BLASTing the genes against the NR database
    and parsing the results form that. We want to blast only the fresh water
    clusters against a modified refseq database."""

    def __init__(self, analysis):
        self.analysis = analysis
        self.clusters = analysis.fresh_clusters

    @property_cached
    def blast_db(self):
        """A blastable database of all refseq + all marine organism genes"""
        db = BLASTdb(refseq_special_db, self.analysis.seq_type)
        if not self.p.all_nin and not self.p.all_pin:
            print "--> STEP 1: Building BLASTable database with all genes..."
            shell_output('gunzip -c %s > %s' % (' '.join(genomes.values()), db))
            assert len(db) == sum(map(len,genomes.values()))
            db.makeblastdb()
            self.timer.print_elapsed()
        return db

