
import sys
from BSB.BSB_Utils.UtilityFunctions import check_python_version, check_package_version
from BSB.BSB_Utils.Launcher import bsb_launch
from BSB.BSB_Utils.Parser import parser


def launch_bsb():

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

    launcher = bsb_launch[arguments.subparser_name]
    launcher(arguments)

def test_func():
    print('yeah')