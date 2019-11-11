import os
from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py
import subprocess

with open("README.md", "r") as fh:
    long_description = fh.read()


def compile_dependency(compilation_command, cwd):
    comp = subprocess.Popen(compilation_command, stdout=subprocess.PIPE,
                            universal_newlines=True, cwd=cwd)
    comp.wait()
    for line in iter(comp.stdout.readline, ''):
        formatted_line = line.strip()
        print(formatted_line)


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
    if not os.path.exists(f'{bowtie2_directory}/bowtie2-align-l'):
        print('Compiling Bowtie2')
        compile_dependency(make_bowtie2_1, bowtie2_directory)
        compile_dependency(make_bowtie2_2, bowtie2_directory)


try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            self.root_is_pure = False
except ImportError:
    bdist_wheel = None


class DevelopCmd(develop):

    def run(self):
        make_external_dependencies()
        super().run()


class BuildCmd(build_py):

    def run(self):
        make_external_dependencies()
        self.data_files = [('BSBolt', 'BSBolt', 'build/lib/BSBolt', []),
                           ('BSBolt.Align', 'BSBolt/Align', 'build/lib/BSBolt/Align', []),
                           ('BSBolt.CallMethylation', 'BSBolt/CallMethylation', 'build/lib/BSBolt/CallMethylation', []),
                           ('BSBolt.External', 'BSBolt/External', 'build/lib/BSBolt/External',
                            ['ART/art_illumina', 'ART/ART_profiler_illumina/art_profiler_illumina',
                             'ART/ART_profiler_illumina/combinedAvg.pl', 'ART/ART_profiler_illumina/empDist.pl',
                             'ART/ART_profiler_illumina/fastqReadAvg.pl', 'ART/ART_profiler_illumina/summation.pl',
                             'ART/Illumina_profiles/Emp100R1.txt', 'ART/Illumina_profiles/Emp100R2.txt',
                             'ART/Illumina_profiles/Emp36R1.txt', 'ART/Illumina_profiles/Emp36R2.txt',
                             'ART/Illumina_profiles/Emp44R1.txt', 'ART/Illumina_profiles/Emp44R2.txt',
                             'ART/Illumina_profiles/Emp50R1.txt', 'ART/Illumina_profiles/Emp50R2.txt',
                             'ART/Illumina_profiles/Emp75R1.txt', 'ART/Illumina_profiles/Emp75R2.txt',
                             'ART/Illumina_profiles/EmpMiSeq250R1.txt', 'ART/Illumina_profiles/EmpMiSeq250R2.txt',
                             'ART/Illumina_profiles/EmpR36R1.txt', 'ART/Illumina_profiles/EmpR36R2.txt',
                             'ART/Illumina_profiles/EmpR44R1.txt', 'ART/Illumina_profiles/EmpR44R2.txt',
                             'ART/Illumina_profiles/EmpR50R1.txt', 'ART/Illumina_profiles/EmpR50R2.txt',
                             'ART/Illumina_profiles/EmpR75R1.txt', 'ART/Illumina_profiles/EmpR75R2.txt',
                             'ART/Illumina_profiles/HiSeq2500L125R1.txt', 'ART/Illumina_profiles/HiSeq2500L125R2.txt',
                             'ART/Illumina_profiles/HiSeq2500L150R1.txt',
                             'ART/Illumina_profiles/HiSeq2500L150R1filter.txt',
                             'ART/Illumina_profiles/HiSeq2500L150R2.txt',
                             'ART/Illumina_profiles/HiSeq2500L150R2filter.txt',
                             'ART/Illumina_profiles/HiSeq2kL100R1.txt', 'ART/Illumina_profiles/HiSeq2kL100R2.txt',
                             'ART/Illumina_profiles/HiSeqXPCRfreeL150R1.txt',
                             'ART/Illumina_profiles/HiSeqXPCRfreeL150R2.txt',
                             'ART/Illumina_profiles/HiSeqXtruSeqL150R1.txt',
                             'ART/Illumina_profiles/HiSeqXtruSeqL150R2.txt', 'ART/Illumina_profiles/MiSeqv3L250R1.txt',
                             'ART/Illumina_profiles/MiSeqv3L250R2.txt', 'ART/Illumina_profiles/MiniSeqTruSeqL50.txt',
                             'ART/Illumina_profiles/NextSeq500v2L75R1.txt',
                             'ART/Illumina_profiles/NextSeq500v2L75R2.txt', 'BT2/bowtie2', 'BT2/bowtie2-build',
                             'BT2/bowtie2-inspect', 'BT2/bowtie2-build-l', 'BT2/bowtie2-build-s',
                             'BT2/bowtie2-inspect-s',
                             'BT2/bowtie2-inspect-l', 'BT2/bowtie2-align-s', 'BT2/bowtie2-align-l']),
                           ('BSBolt.Impute', 'BSBolt/Impute', 'build/lib/BSBolt/Impute', []),
                           ('BSBolt.Impute.Imputation', 'BSBolt/Impute/Imputation',
                            'build/lib/BSBolt/Impute/Imputation', []),
                           ('BSBolt.Impute.Impute_Utils', 'BSBolt/Impute/Impute_Utils',
                            'build/lib/BSBolt/Impute/Impute_Utils', []), (
                           'BSBolt.Impute.Validation', 'BSBolt/Impute/Validation', 'build/lib/BSBolt/Impute/Validation',
                           []),
                           ('BSBolt.Index', 'BSBolt/Index', 'build/lib/BSBolt/Index', []),
                           ('BSBolt.Matrix', 'BSBolt/Matrix', 'build/lib/BSBolt/Matrix', []),
                           ('BSBolt.Simulate', 'BSBolt/Simulate', 'build/lib/BSBolt/Simulate', []),
                           ('BSBolt.Utils', 'BSBolt/Utils', 'build/lib/BSBolt/Utils', []),
                           ('BSBolt.Variant', 'BSBolt/Variant', 'build/lib/BSBolt/Variant', [])]
        super().run()


setup(name='BSBolt',
      version='0.1.2',
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
      classifiers=['Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8'],
      platforms=["Linux", "Mac OS-X", "Unix"],
      requires=['pysam', 'numpy', 'tqdm'],
      install_requires=['pysam==0.15.2', 'numpy>=1.16.3', 'tqdm>=4.31.1'],
      entry_points={'console_scripts': ['BSBolt = BSBolt.__main__:launch_bsb']},
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
      cmdclass={'develop': DevelopCmd,
                'build_py': BuildCmd,
                'bdist_wheel': bdist_wheel},
      zip_safe=False
      )
