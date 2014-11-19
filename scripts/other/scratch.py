from ld12.analysis import Analysis; out_dir = "/Users/sinclair/proj/b2013274/results"; a = Analysis(out_dir=out_dir)

with open(a.p.user_dir + 'best_clusters.txt', 'w') as handle: handle.writelines('\n'.join([c.name for c in a.best_clusters]))

with open(a.p.user_dir + 'fresh_clusters.txt', 'w') as handle: handle.writelines('\n'.join([c.name for c in a.fresh_clusters]))

from tqdm import tqdm
for c in tqdm(a.fresh_clusters): c.fasta