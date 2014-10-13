"""How are we going to build the ribosomal tree"""

###############################################################################
from ld12 import genes, genomes

ribo_genes = [g for g in genes.values() if g.ribo_group]
print len(ribo_genes)

ribo_groups = set([g.ribo_group for g in ribo_genes])
print len(ribo_groups)

###############################################################################
from collections import defaultdict
import pandas

result = defaultdict(lambda: defaultdict(int))
for group in ribo_groups:
    current_genes = [g for g in ribo_genes if g.ribo_group==group]
    for gene in current_genes:
        result[group][gene.genome.name] += 1

result = pandas.DataFrame(result, dtype=int)
result = result.fillna(0)
result = result.astype(int)

###############################################################################
from ld12 import repos_dir
result.to_csv(repos_dir + "scripts/ribo_table.tsv", sep='\t', encoding='utf-8')