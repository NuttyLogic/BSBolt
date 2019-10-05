#! /usr/bin/env python3

import sys
from BSB.BSB_Utils.UtilityFunctions import check_python_version, get_external_paths, check_package_version
from BSB.BSB_Utils.Launcher import bsb_launch
from BSB.BSB_Utils.Parser import parser

bt2_path, art_path = get_external_paths()
check_python_version()
if not check_package_version():
    print('Please Update Package Requirements')
    print('cd BSBolt')
    print('pip3 install -r requirements.txt')
    sys.exit()

arguments = parser.parse_args()

if len(sys.argv[1:]) == 0:
    parser.print_help()
    # parser.print_usage() # for just the usage line
    parser.exit()


if __name__ == "__main__":
    launcher = bsb_launch[arguments.subparser_name]
    launcher(arguments)
