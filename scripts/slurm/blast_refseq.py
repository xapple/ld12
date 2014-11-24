#!/usr/bin/env python
#SBATCH -D /home/lucass/LD12/
#SBATCH -J ld12_blastp_refseq
#SBATCH -o /home/lucass/LD12/blastp_refseq.out
#SBATCH -A b2011105
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
now = now.strftime("%Y-%m-%d %Hh%Mm%Ss %Z%z")
print "SLURM: start at {0} on {1}".format(now, platform.node())

# Do it #
a = Analysis(out_dir=out_dir)
a.duplications.assign_best_hits()
a.duplications.assign_taxonomy()
a.duplications.save_duplications_stats()

# End #
now = datetime.datetime.now(dateutil.tz.tzlocal())
now = now.strftime("%Y-%m-%d %Hh%Mm%Ss %Z%z")
print "SLURM: end at {0}".format(now)