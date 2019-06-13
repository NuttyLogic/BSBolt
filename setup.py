from distutils.core import setup

setup(name='BiSulfiteBolt',
      version='0.2.0',
      description='Bisulfite Sequencing Processing Platform',
      author='Colin Farrell',
      author_email='colinpfarrell@gmail.com',
      license='GPLv3',
      packages=['BSBolt'],
      requires=['pysam', 'numpy', 'tqdm'])
