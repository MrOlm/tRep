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
    stb = args.get('scaffold2bin', None)

    # Make Tdb (gene_level taxonomy)
    Tdb = tRep.controller.convert_b6_to_Tdb(args, save=False)
    Tdb.to_csv(os.path.join(out_base) + '_fullGeneTaxonomy.tsv', \
        index=False, sep='\t')

    # Make genome level taxonomy
    if stb is None:
        pass
    else:
        if stb == 'ALL':
            Tdb['bin'] = 'genome'
        else:
            Tdb = tRep.add_bin_to_tdb(Tdb, stb)
        gdb = tRep.gen_taxonomy_table(Tdb, on='bin')
        gdb.to_csv(os.path.join(out_base) + '_fullGenomeTaxonomy.tsv', \
            index=False, sep='\t')

    # Make scaffold level taxonomy
    if not skip_scaffs:
        Tdb['scaffold'] = ['_'.join(x.split('|')[0].split('_')[:-1]) for x in Tdb['querry']]
        try:
            sdb = tRep.gen_taxonomy_table(Tdb, on='scaffold')
            sdb.to_csv(os.path.join(out_base) + '_fullScaffoldTaxonomy.tsv', \
                index=False, sep='\t')
        except:
            print('unable to parse scaffold information- skipping')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(\
        description=\
'''
Generate taxonomy report from a b6+ file
''',\
        formatter_class=argparse.RawTextHelpFormatter)

    InpArgs = parser.add_argument_group('INPUT ARGUMENTS')
    InpArgs.add_argument('-b', '--b6_loc', help='location of b6+ file', required=True)
    InpArgs.add_argument('-stb', '--scaffold2bin',  help=\
        'scaffold to bin file for generating per-genome taxonomy. Pass in the word "ALL" to consider all genes the same genome')
    InpArgs.add_argument('-a', '--amino_acids',  help=\
        'prodigal output used to make the b6 file; allows identification of fully unknown genes')

    OutArgs = parser.add_argument_group('OUTPUT ARGUMENTS')
    OutArgs.add_argument('-o', '--out_loc',  help=\
        'output basename', required=True)

    OptArgs = parser.add_argument_group('OPTIONAL ARGUMENTS')
    OptArgs.add_argument('--SkipScaffolds',  help=\
        'skip generating per-scaffold taxonomy',action='store_true', default=False)
    OptArgs.add_argument('--tax_type',  help=\
        "If using the translated b6, do you want the species or group taxID?",
        action='store', default='species', choices=['species', 'group'])



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
