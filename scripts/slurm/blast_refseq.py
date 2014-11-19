#!/usr/bin/env python
#SBATCH -D /home/lucass/LD12/
#SBATCH -J ld12_tblastn_refseq
#SBATCH -o /home/lucass/LD12/ld12_tblastn_refseq.out
#SBATCH -A b2011138
#SBATCH -t 7-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user lucas.sinclair@me.com
#SBATCH --mail-type=END
#SBATCH -d singleton

# Modules #
import dateutil, datetime, platform, socket, os
from ld12.analysis import Analysis

# Constants #
host = socket.gethostname()
home = os.environ['HOME'] + '/'

# Output path #
if host.startswith('milou'): out_dir = home + "/proj/b2013274/results"
else:                        out_dir = home + "/LD12/results"

# Message #
now = datetime.datetime.now(dateutil.tz.tzlocal())
now = now.strftime("%Y-%m-%d %H:%M:%S %Z%z")
print "SLURM: start at {0} on {1}".format(now, platform.node())

# Do it #
a = Analysis(out_dir=out_dir)
print a.duplications.search_results

# End #
now = datetime.datetime.now(dateutil.tz.tzlocal())
now = now.strftime("%Y-%m-%d %H:%M:%S %Z%z")
print "SLURM: end at {0}".format(now)