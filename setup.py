from distutils.core import setup

setup(
      name             = 'ld12',
      version          = '1.0.1',
      description      = 'Analysis for LD12 genomes and the freshwater against marine comparison paper.',
      long_description = open('README.md').read(),
      license          = 'MIT',
      url              = 'https://github.com/xapple/ld12',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['ld12'],
      scripts          = ['ld12/ld12'],
      requires         = ['plumbing', 'fasta', 'seqsearch', 'biopython', 'sh', 'pandas', 'shell_command'],
)