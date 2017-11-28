#!/usr/bin/env python

import argparse
import sys
import os

import tRep
import tRep.controller

__version__ = tRep.__version__

def main(**args):
    out_base = args.get('out_loc')

    Tdb = tRep.controller.convert_b6_to_Tdb(args, save=False)
    Tdb.to_csv(os.path.join(out_base, 'Tdb'))

    try:
        Tdb['genome'] = ''

if __name__ == "__main__":
    parser = argparse.ArgumentParser(\
        description='''
Generate taxonomy report from b6+ file
                    ''',\
        formatter_class=argparse.RawTextHelpFormatter)

    InpArgs = parser.add_argument_group('INPUT ARGUMENTS')
    InpArgs.add_argument('-b', '--b6_loc', help='location of b6+ file')

    OutArgs = parser.add_argument_group('OUTPUT ARGUMENTS')
    OutArgs.add_argument('-o', '--out_loc',  help=\
        'output basename')

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = sys.argv[1:]
    if (len(args) == 0):
        print('Run with -h for help')
        sys.exit(0)

    args = parser.parse_args()
    main(**vars(args))
