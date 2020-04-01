#! /usr/bin/env python3

import sys
from BSBolt.Utils.UtilityFunctions import check_python_version, check_package_version
from BSBolt.Utils.Launcher import bsb_launch
from BSBolt.Utils.Parser import parser


def launch_bsb():
    check_python_version()
    if not check_package_version():
        print('Please Update Packages')
        sys.exit()

    arguments = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()

    launcher = bsb_launch[arguments.subparser_name]
    launcher(arguments)


if __name__ == "__main__":
    launch_bsb()
