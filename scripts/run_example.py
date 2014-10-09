import socket, os
host = socket.gethostname()
home = os.environ['HOME'] + '/'

if host.startswith('milou'): out_dir = home + "/proj/b2013274/results"
else:                        out_dir = home + "/proj/b2013274/results"

from ld12.analysis import Analysis
a = Analysis(out_dir=out_dir)
a.search_results