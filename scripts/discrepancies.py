"""There are some discrepancies between the fna and faa files. What are they"""

import numpy
from fasta import FASTA

names = ['2236446587', '2236446227', '2236446226', '2236446115', '2236446114', '2236445762', '2236446248', '2236446050', '2236446154', '2236446074', '2236445702', '2236446303', '2236446345', '2236446625', '2236446051', '2236446049', '2236446581', '2236445824', '2236446355', '2236446351', '2236446370', '2236446339', '2236446062', '2236445915', '2236445916', '2236445669', '2236446147', '2236446048', '2236446182']
path = "/proj/b2013274/mcl/2236347014.genes.fna"
fasta = FASTA(path)

print "all:", numpy.mean(fasta.lengths)
print "missing:", numpy.mean([len(fasta[name]) for name in names])
