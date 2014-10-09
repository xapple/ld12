import socket, os
host = socket.gethostname()
home = os.environ['HOME'] + '/'

if host.startswith('milou'): out_dir = home + "/proj/b2013274/results"
else:                        out_dir = home + "/proj/b2013274/results"

from ld12.analysis import Analysis
a = Analysis(out_dir=out_dir)
a.save_count_table()
for c in a.best_clusters[0:2]: print c.phylogeny