from distutils.core import setup

setup(
      name             = 'ld12',
      version          = '0.0.1',
      description      = 'Analysis for LD12 genomes.',
      long_description = open('README.md').read(),
      license          = 'MIT',
      url              = 'https://bitbucket.org/xapple/ld12',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['ld12'],
      requires         = ['plumbing', 'fasta', 'seqsearch', 'biopython', 'sh', 'pandas', 'shell_command'],
)