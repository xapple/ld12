from ld12.analysis import Analysis; out_dir = "/Users/sinclair/proj/b2013274/results"; a = Analysis(out_dir=out_dir)

with open(a.p.user_dir + 'best_clusters.txt', 'w') as handle: handle.writelines('\n'.join([c.name for c in a.best_clusters]))

with open(a.p.user_dir + 'fresh_clusters.txt', 'w') as handle: handle.writelines('\n'.join([c.name for c in a.fresh_clusters]))

from tqdm import tqdm
for c in tqdm(a.fresh_clusters): c.fasta

# Count the top alphaproteobcteria (besides SAR11) #
a = Analysis()
a.duplications.assign_hits()
from tqdm import tqdm
from collections import Counter
counts_class = Counter()
counts_full = Counter()

def process(genes):
    for g in tqdm(genes):
        tax = g.hits[-1]['taxonomy'].split(';')
        tax += ["unassigned", "unassigned", "unassigned"]
        if   tax[0].strip() != 'Bacteria':             continue
        elif tax[1].strip() != 'Proteobacteria':       continue
        elif tax[2].strip() != 'Alphaproteobacteria':  continue
        else:
            counts_class.update([tax[3].strip()])
            counts_full.update([tuple(tax)])

process(a.duplications.best_is_other)