from distutils.version import LooseVersion
import os
import sys
import pysam


def reverse_complement(sequence):
    """
    Params:

    * *sequence (str)*: DNA sequence, non ATGC nucleotide will be returned unaltered

    Returns:

    * *reversed_string.translate(_rc_trans) (str)*: reverse complement of input sequence
    """
    # reverse string
    reversed_string = sequence[::-1]
    # replace base with complement
    _rc_trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return reversed_string.translate(_rc_trans)


def retrieve_iupac(nucleotide):
    """
    Params:

    * *nucleotide (str)*: single character

    Returns:

    * *iupac_tuple (tuple)*: tuple of strings with possible bases
    """
    iupac_key = {'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'), 'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                 'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'), 'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G'),
                 'N': ('A', 'C', 'G', 'T')}
    try:
        iupac_tuple = iupac_key[nucleotide]
    except KeyError:
        iupac_tuple = tuple(nucleotide)
    return iupac_tuple


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
        if LooseVersion(installed_version) < LooseVersion(required_version):
            print(f'{package} {installed_version} < {package} {required_version}, please update {package}')
            print(f'pip3 install {package} --upgrade')
            all_packages = False
    return all_packages


def check_python_version():
    if sys.version_info < (3, 6):
        print('Python must be >= 3.6.0')
        raise OSError


def get_external_paths():
    """Get paths of dependencies. Print warning if setup.py not run and dependencies not compiled.

    Returns:
    * *bwa (str)*: path to bwa executable
    * *wgsim (str)*: path to wgsim executable
    """
    utility_directory = os.path.dirname(os.path.realpath(__file__))
    external_directory = '/'.join(utility_directory.split('/')[:-1]) + '/External/'
    bwa = f'{external_directory}BWA/bwa-mem2'
    wgsim = f'{external_directory}WGSIM/wgsim'
    if not os.path.exists(bwa) or not os.path.exists(wgsim):
        print(f'Must compile external dependencies\n python3 setup.py build')
    return bwa, wgsim


def propagate_error(error):
    raise error


def sort_bam(bam_output=None, bam_input=None):
    """ Sort bam file

    Params:
    * *bam_output (str)*: output path for sorted bam file
    * *bam_input (str)*: input bam file"""
    pysam.sort('-o', bam_output, bam_input)


def index_bam(bam_input=None):
    """
    Index bam file

    Params:
    * *bam_input (str)*: input bam files
    """
    pysam.index(bam_input)
