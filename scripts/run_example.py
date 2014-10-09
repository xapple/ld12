# Modules #
from tqdm import tqdm
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
for c in tqdm(a.best_clusters): c.tree