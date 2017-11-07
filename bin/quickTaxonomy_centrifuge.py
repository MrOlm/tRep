#!/usr/bin/env python

"""
tRep
"""
import argparse
import sys

import tRep

__version__ = tRep.__version__

if __name__ == "__main__":
    """ This is executed when run from the command line
    """
    parser = argparse.ArgumentParser(description="Quickly determine the taxonomy of a genome")

    InpArgs = parser.add_argument_group('INPUT ARGUMENTS')
    # Parsing from raw centrifuge
    InpArgs.add_argument('-ch', "--cent_hit",\
        help='raw output of centrifuge run (hits.tsv)')
    InpArgs.add_argument('-cr', "--cent_report",\
        help='raw output of centrifuge run (report.tsv)')
    # Parsing from Tdb (from dRep)
    # InpArgs.add_argument('-tdb',
    #     help='Tdb.csv (from dRep run))
    # Parsing from .fasta file
    InpArgs.add_argument('-f', "--fasta",
        help='location of fasta file(s) to detemine the taxonomy of')
    InpArgs.add_argument("--cent_index",help='path to centrifuge index (for example, ' + \
                    "/home/mattolm/download/centrifuge/indices/b+h+v", default= \
                    "/home/mattolm/download/centrifuge/indices/b+h+v")

    OutArgs = parser.add_argument_group('OUTPUT ARGUMENTS')
    OutArgs.add_argument('-m', '--method', default='percent', help=\
        'Method to determine taxonomy. percent = the lowest taxonomic id with ' + \
        ' at least -p hits. max = the taxonomic level with the most hits (lowest' + \
        ' is genus')
    OutArgs.add_argument('-per', '--percent', default='50', help=\
        'Minimum percent of genes for the percent method')
    OutArgs.add_argument('--full_dump', default=None, help=\
        'Output FULL taxonomy table to the file specified here (for testing purposes)')

    FilArgs = parser.add_argument_group('FILTERING ARGUMENTS')
    FilArgs.add_argument('--min_score', default='250', help=\
        'Minimum score for centrifuge hit')
    FilArgs.add_argument('--min_diff', default='0', help=\
        'Minimum score difference between first and second hit')

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))
    parser.add_argument('-p', '--processes', help='threads for centrifuge run', default='6')

    args = sys.argv[1:]
    if (len(args) == 0):
        print('Run with -h for help')
        sys.exit(0)

    args = parser.parse_args()
    tRep.main(**vars(args))
