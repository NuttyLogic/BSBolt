import os
from setuptools import Extension, setup
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py
import subprocess

with open("README.md", "r") as fh:
    long_description = fh.read()

"""
External dependencies compiled for linux using quay.io/pypa/manylinux2010_x86_64 image:
/opt/python/cpXX-cpXX/bin/python setup.py bdist_wheel
"""

def compile_dependency(compilation_command, cwd):
    comp = subprocess.Popen(compilation_command, stdout=subprocess.PIPE,
                            universal_newlines=True, cwd=cwd)
    comp.wait()
    for line in iter(comp.stdout.readline, ''):
        formatted_line = line.strip()
        print(formatted_line)


def make_external_dependencies():
    working_directory = os.path.dirname(os.path.realpath(__file__))
    bwa_directory = f'{working_directory}/BSBolt/External/BWA'
    wgsim_directory = f'{working_directory}/BSBolt/External/WGSIM'
    if not os.path.exists(f'{wgsim_directory}/wgsim'):
        print('Compiling wgsim')
        compile_dependency(['make', '-portable'], wgsim_directory)
    if not os.path.exists(f'{bwa_directory}/bwa-mem2'):
        print('Compiling bwa-mem2')
        compile_dependency(['make'], bwa_directory)


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
                            ['WGSIM/wgsim', 'BWA/bwa-mem2']),
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
      version='0.3.0',
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
      install_requires=['pysam>=0.15.4', 'numpy>=1.16.3', 'tqdm>=4.31.1'],
      entry_points={'console_scripts': ['BSBolt = BSBolt.__main__:launch_bsb']},
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
      cmdclass={'develop': DevelopCmd,
                'build_py': BuildCmd,
                'bdist_wheel': bdist_wheel},
      zip_safe=False
      )
