#!/usr/bin/env python

import argparse
import sys
import os

import tRep
import tRep.controller

__version__ = tRep.__version__

def main(**args):
    out_base = args.get('out_loc')
    skip_scaffs = args.get('SkipScaffolds', False)

    # Make Tdb
    Tdb = tRep.controller.convert_b6_to_Tdb(args, save=False)
    Tdb.to_csv(os.path.join(out_base) + '_fullGeneTaxonomy.tsv', \
        index=False, sep='\t')

    # Make genome level taxonomy
    Tdb['genome'] = 'genome'
    gdb = tRep.gen_taxonomy_table(Tdb, on='genome')
    gdb.to_csv(os.path.join(out_base) + '_fullGenomeTaxonomy.csv', \
        index=False, sep='\t')

    # Make scaffold level taxonomy
    if not skip_scaffs:
        try:
            Tdb['scaffold'] = ['_'.join(x.split('_')[:-1]) for x in Tdb['querry']]
            sdb = tRep.gen_taxonomy_table(Tdb, on='scaffold')
            sdb.to_csv(os.path.join(out_base) + '_fullScaffoldTaxonomy.tsv', \
                index=False, sep='\t')
        except:
            print('unable to parse scaffold information- skipping')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(\
        description=\
'''
Generate taxonomy report from a b6+ file resulting from a single genome
''',\
        formatter_class=argparse.RawTextHelpFormatter)

    InpArgs = parser.add_argument_group('INPUT ARGUMENTS')
    InpArgs.add_argument('-b', '--b6_loc', help='location of b6+ file')

    OutArgs = parser.add_argument_group('OUTPUT ARGUMENTS')
    OutArgs.add_argument('-o', '--out_loc',  help=\
        'output basename')

    OptArgs = parser.add_argument_group('OPTIONALL ARGUMENTS')
    OptArgs.add_argument('--SkipScaffolds',  help=\
        'skip generating per-scaffold taxonomy',action='store_true')
    OptArgs.add_argument('-stb', '--scaffold2bin',  help=\
        'scaffold to bin file for generating per-genome taxonomy')

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
