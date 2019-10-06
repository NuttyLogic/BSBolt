import os
from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
import subprocess

with open("README.md", "r") as fh:
    long_description = fh.read()
working_directory = os.path.dirname(os.path.realpath(__file__))
bowtie2_directory = f'{working_directory}/BSBolt/External/BT2'
art_directory = f'{working_directory}/BSBolt/External/ART'


def compile_dependency(compilation_command, cwd):
    print(' '.join(compilation_command))
    comp = subprocess.Popen(compilation_command, stdout=subprocess.PIPE,
                            universal_newlines=True, cwd=cwd)
    for line in iter(comp.stdout.readline, ''):
        print(line)


def make_external_dependencies():
    make_bowtie2_1 = ['make', 'static-libs']
    make_bowtie2_2 = ['make', 'STATIC_BUILD=1']
    art_src = 'art_illumina_src/'
    make_art = ['g++', '-static', '-g', '-O2', '-o', 'art_illumina', f'{art_src}art_illumina.o',
                f'{art_src}art_qual_scale.o', f'{art_src}empdist.o', f'{art_src}readSeqFile.o',
                f'{art_src}seqRead.o', f'{art_src}samRead.o', '-lgsl', '-lgslcblas', '-lm']
    if not os.path.exists(f'{art_directory}/art_illumina'):
        print('Compiling ART')
        compile_dependency(make_art, art_directory)
    if not os.path.exists(f'{bowtie2_directory}/bowtie2-align-l'):
        print('Compiling Bowtie2')
        compile_dependency(make_bowtie2_1, bowtie2_directory)
        compile_dependency(make_bowtie2_2, bowtie2_directory)


class InstallCmd(install):

    def run(self):
        make_external_dependencies()
        install.run(self)


class DevelopCmd(develop):

    def run(self):
        make_external_dependencies()
        develop.run(self)


make_external_dependencies()

setup(name='BSBolt',
      version='0.1.1',
      description='Bisulfite Sequencing Processing Platform',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/BSBolt',
      project_urls={'Documentation': 'https://bsbolt.readthedocs.io/en/latest/'},
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='GPLv3',
      packages=['BSBolt',
                'BSBolt.Align',
                'BSBolt.CallMethylation',
                'BSBolt.External',
                'BSBolt.Impute',
                'BSBolt.Impute.Imputation',
                'BSBolt.Impute.Impute_Utils',
                'BSBolt.Impute.Validation',
                'BSBolt.Index',
                'BSBolt.Matrix',
                'BSBolt.Simulate',
                'BSBolt.Utils',
                'BSBolt.Variant'],
      requires=['pysam', 'numpy', 'tqdm'],
      install_requires=['pysam>=0.15.2', 'numpy>=1.16.3', 'tqdm>=4.31.1'],
      entry_points={'console_scripts': ['BSBolt = BSBolt.Launch:launch_bsb']},
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True
      )
