outdir = "/proj/b2013274/results"
outdir = "/glob/lucass/other/alex_protein/"

from ld12.analysis import Analysis
a = Analysis(outdir=outdir)
a.run()


###############################################################################
output_directory = home + "glob/lucass/other/alex_protein/"
extension = '.faa'

# Real input #
files = "/proj/b2013274/mcl/*" + extension
genomes = [Genome(path) for path in glob.glob(files)]
analysis = Analysis(genomes, output_directory)

# Test input #
test_files = output_directory + "/test/*" + extension
test_genomes = [Genome(path) for path in glob.glob(test_files)]
test_analysis = Analysis(test_genomes, output_directory + "/test/")

# Main program #
def run(): Phylo.draw_ascii(analysis.master_cluster.phylogeny)
