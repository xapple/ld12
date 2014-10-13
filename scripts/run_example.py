"""You can use the command line script to run the pipeline, or you can use this example to type commands into an ipython prompt yourself"""

# Modules #
import socket, os
from ld12.analysis import Analysis

# Constants #
host = socket.gethostname()
home = os.environ['HOME'] + '/'

# Output path #
if host.startswith('milou'): out_dir = home + "/proj/b2013274/results"
else:                        out_dir = home + "/proj/b2013274/results"

# Do it #
a = Analysis(out_dir=out_dir)
a.make_trees()