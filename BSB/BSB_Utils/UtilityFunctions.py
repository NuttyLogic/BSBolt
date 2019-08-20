from distutils.version import LooseVersion
import os
import platform
import subprocess
import sys
from packaging import version


def reverse_complement(sequence):
    """
    Arguments:
        sequence (str): DNA sequence, can have non ATGC nucleotide will remain untouched
    Returns:
       reversed_string.translate(_rc_trans) (str): reverse complement of input sequence
    """
    # reverse string
    reversed_string = sequence[::-1]
    # replace base with complement
    _rc_trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return reversed_string.translate(_rc_trans)


def retrieve_iupac(nucleotide):
    """
    Arguments:
        nucleotide (str): single character
    Returns:
        iupac_tuple (tuple): tuple of strings with possible bases
    """
    iupac_key = {'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'), 'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                 'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'), 'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G'),
                 'N': ('A', 'C', 'G', 'T')}
    try:
        iupac_tuple = iupac_key[nucleotide]
    except KeyError:
        iupac_tuple = tuple(nucleotide)
    return iupac_tuple


def check_bowtie2_path(bowtie2_path='bowtie2'):
    """Simple function to check bowti2 path. If path is valid and version >= 2.3.4.2 function exits normally. If path
    not valid raise FileNotFoundError. If version < 2.3.4.2 raise RuntimeWarning.
    Keyword Arguments:
        bowtie2_path (str): path to bowtie executable, default = bowtie2, default assumes bowtie2 is in system path
    """
    # external command
    bowtie2_version_command = [bowtie2_path, '--version']
    try:
        bowtie2_check = subprocess.Popen(bowtie2_version_command, stdout=subprocess.PIPE, universal_newlines=True)
    except FileNotFoundError:
        raise FileNotFoundError('Bowtie2 not in system path or Bowtie2 Executable Path Incorrect')
    else:
        # get first line of stdout
        version_line = next(iter(bowtie2_check.stdout.readline, ''))
        bowtie2_version = version_line.replace('\n', '').split(' ')[-1]
        if LooseVersion(bowtie2_version) < LooseVersion('2.2.9'):
            raise RuntimeWarning('BSeeker-R Performance not evaluated on Bowtie2 Versions < 2.2.9')


def import_package_check(package_name):
    try:
        package = __import__(package_name)
    except ModuleNotFoundError:
        print(f'{package_name} not found, please install {package_name}')
        print(f'pip3 install {package_name} --user')
        sys.exit()
    else:
        return package.__version__


def check_package_version():
    """Check third party packages to make sure version is greater than requirement"""
    packages = dict(pysam='0.15.2', tqdm='4.31.1', numpy='1.16.3')
    all_packages = True
    for package, required_version in packages.items():
        installed_version = import_package_check(package)
        if version.parse(installed_version) < version.parse(required_version):
            print(f'{package} {installed_version} < {package} {required_version} requirement, please update')
            print(f'pip3 install {package} --upgrade')
            all_packages = False
    return all_packages


def check_python_version():
    if sys.version_info < (3, 6):
        print('Python must be >= 3.6.0')
        raise OSError


def get_external_paths():
    utility_directory = os.path.dirname(os.path.realpath(__file__))
    external_directory = '/'.join(utility_directory.split('/')[:-2]) + '/External/'
    bt2 = None
    art = None
    if platform.system() == 'Linux':
        bt2 = f'{external_directory}BT2/bowtie2-2.3.4.3-linux-x86_64/bowtie2'
        art = f'{external_directory}ART/ART_linux/art_illumina'
        return bt2, art
    elif platform.system() == 'Darwin':
        bt2 = f'{external_directory}BT2/bowtie2-2.3.4.3-macos-x86_64/bowtie2'
        art = f'{external_directory}ART/ART_mac/art_illumina'
        return bt2, art
    if not bt2 and not art:
        print(f'Unsupported platform version; {platform.system()}. Please proved bowtie2 paths and ART paths. '
              f'BSBolts functionality may be affected.')
