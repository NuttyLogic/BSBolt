import os
from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_ext import build_ext
import subprocess
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()


def compile_dependency(compilation_command, cwd):
    comp = subprocess.Popen(compilation_command, stdout=subprocess.PIPE,
                            universal_newlines=True, cwd=cwd)
    comp.wait()
    for line in iter(comp.stdout.readline, ''):
        print(line)


def make_external_dependencies():
    working_directory = os.path.dirname(os.path.realpath(__file__))
    bowtie2_directory = f'{working_directory}/BSBolt/External/BT2'
    art_directory = f'{working_directory}/BSBolt/External/ART'
    make_bowtie2_1 = ['make', 'static-libs']
    make_bowtie2_2 = ['make', 'STATIC_BUILD=1']
    make_art = ['sh', 'art_comp.sh']
    if not os.path.exists(f'{art_directory}/art_illumina'):
        print('Compiling ART')
        compile_dependency(make_art, art_directory)
        print(os.path.exists(f'{art_directory}/art_illumina'))
        print(glob.glob(f'{art_directory}/*'))
    if not os.path.exists(f'{bowtie2_directory}/bowtie2-align-l'):
        print('Compiling Bowtie2')
        compile_dependency(make_bowtie2_1, bowtie2_directory)
        compile_dependency(make_bowtie2_2, bowtie2_directory)


class InstallCmd(install):

    def run(self):
        print('YESSLKSDLKJFJKL:SKL:DJFKL:J')
        make_external_dependencies()
        super().run()


class DevelopCmd(develop):

    def run(self):
        make_external_dependencies()
        super().run()


class BuildCmd(build_ext):

    def run(self):
        make_external_dependencies()
        super().run()



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
      include_package_data=True,
      cmdclass={'install': InstallCmd,
                'develop': DevelopCmd,
                'build_ext': BuildCmd},
      zip_safe=False
      )
