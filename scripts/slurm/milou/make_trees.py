#!/usr/bin/env python
#SBATCH -D /home/lucass/LD12/
#SBATCH -J ld12_make_trees
#SBATCH -o /home/lucass/LD12/make_trees.out
#SBATCH -A b2011138
#SBATCH -t 7-00:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user lucas.sinclair@me.com
#SBATCH --mail-type=END
#SBATCH -d singleton

# Modules #
import dateutil, datetime, platform, os
from ld12.analysis import Analysis

# Constants #
home    = os.environ['HOME'] + '/'
out_dir = home + "/LD12/results"

# Message #
now = datetime.datetime.now(dateutil.tz.tzlocal())
now = now.strftime("%Y-%m-%d %H:%M:%S %Z%z")
print "SLURM: start at {0} on {1}".format(now, platform.node())

# Do it #
a = Analysis(out_dir=out_dir)

# Long job #
#a.make_trees()
assert all([c.p.bestTree.exists for c in a.best_clusters])

# Stats #
a.comparison.save_uncollapsible_stats()
a.comparison.save_split_three_a_b()
a.comparison.save_split_conserved()
a.comparison.save_mismatching_stats()

# End #
now = datetime.datetime.now(dateutil.tz.tzlocal())
now = now.strftime("%Y-%m-%d %H:%M:%S %Z%z")
print "SLURM: end at {0}".format(now)