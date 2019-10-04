from setuptools import setup, find_packages, find_namespace_packages

setup(name='BSBolt',
      version='0.1.0',
      description='Bisulfite Sequencing Processing Platform',
      url='https://github.com/NuttyLogic/BSBolt',
      project_urls={'Documentation': 'https://bsbolt.readthedocs.io/en/latest/'},
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='GPLv3',
      packages=find_namespace_packages(),
      requires=['pysam', 'numpy', 'tqdm'],
      entry_points={'console_scripts': ['BSBolt = BSB.RunBSB:test_func']},
      test_suite='Test.test')
