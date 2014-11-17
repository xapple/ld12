#!/usr/bin/env python
#SBATCH -D /home/lucass/LD12/
#SBATCH -J ld12_make_trees
#SBATCH -o /home/lucass/LD12/run.out
#SBATCH -A b2011138
#SBATCH -t 7-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user lucas.sinclair@me.com
#SBATCH --mail-type=END
#SBATCH -d singleton

# Modules #
import time, platform, socket, os
from ld12.analysis import Analysis

# Constants #
host = socket.gethostname()
home = os.environ['HOME'] + '/'

# Output path #
if host.startswith('milou'): out_dir = home + "/proj/b2013274/results"
else:                        out_dir = home + "/LD12/results"

# Do it #
print "SLURM: start at {0} on {1}".format(time.asctime(), platform.node())
a = Analysis(out_dir=out_dir)
a.make_trees()
print "SLURM: end at {0}".format(time.asctime())