"""We explore the client given inputs, check for problems,
then format them and store them in the repository as immutable text files"""

import inspect, os, glob, pandas
from fasta import FASTA

current_script = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(current_script)) + '/'
genomes_dir = current_dir + '../ld12/data/genomes/'

input_dir = "/proj/b2013274/mcl/"
faa_paths = sorted(glob.glob(input_dir + '*.faa'))
fna_paths = sorted(glob.glob(input_dir + '*.fna'))
faas = [FASTA(faa) for faa in faa_paths if '647533246' not in faa]
fnas = [FASTA(fna) for fna in fna_paths if '647533246' not in fna]

faas_nums = [int(g.short_prefix) for g in faas]
fnas_nums = [int(g.short_prefix) for g in fnas]
metadata = pandas.io.parsers.read_csv(current_dir + '../ld12/data/metadata.tsv', sep='\t', index_col=0, encoding='utf-8')
meta_nums = list(metadata.index)

print set(faas_nums) ^ set(fnas_nums)
print set(faas_nums) ^ set(meta_nums)

def strip(seq):
    seq = seq.id
    seq = seq.split('|')[0]
    seq = seq.split('[')[0]
    return seq

for faa,fna in zip(faas, fnas):
    faas_genes = [strip(seq) for seq in faa]
    fnas_genes = [strip(seq) for seq in fna]
    print faa, len(set(fnas_genes) ^ set(faas_genes))

fnas_genes = [strip(seq) for fna in fnas for seq in fna]
print len(fnas_genes), len(set(fnas_genes))

for genome in faas:
    out_path = genomes_dir + genome.short_prefix + '.fasta'
    out_fasta = FASTA(out_path)
    out_fasta.create()
    for seq in genome: out_fasta.add_str(str(seq.seq), strip(seq))
    out_fasta.close()
    out_fasta.gzip_to()
    out_fasta.remove()

def lines():
    for genome in faas:
        for gene in genome:
            name = strip(gene)
            yield name + '\t' + gene.description[len(name):].rstrip(' |') + '\n'

annotations_path = current_dir + '../ld12/data/annotations.tsv'
with open(annotations_path, 'w') as handle: handle.writelines(lines())