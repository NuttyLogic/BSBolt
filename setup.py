from distutils.core import setup

setup(name='BiSulfiteBolt',
      version='0.0.1',
      description='Bisulfite Sequencing Processing Platform',
      author='Colin Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      packages=['BSBolt'],
      requires=['pysam', 'numpy', 'tqdm'])
