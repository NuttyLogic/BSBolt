from setuptools import setup, find_packages, find_namespace_packages

setup(name='BSBolt',
      version='0.1.1',
      description='Bisulfite Sequencing Processing Platform',
      url='https://github.com/NuttyLogic/BSBolt',
      project_urls={'Documentation': 'https://bsbolt.readthedocs.io/en/latest/'},
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='GPLv3',
      packages=['BSB'],
      dependencies=['pysam >= 0.15.2', 'tqdm >= 4.31.1', 'numpy >= 1.16.3'],
      requires=['pysam', 'numpy', 'tqdm'],
      entry_points={'console_scripts': ['BSBolt = BSB.BSBRun:launch_bsb']},
      python_requires='>=3.6',
      test_suite='Test.test',
      include_package_data=True)
