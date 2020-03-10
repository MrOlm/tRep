#!/usr/bin/env python

'''
Simple script to get functional taxonomy from a b6 file and a translation table
'''

import os
import sys
import time
import gzip
import argparse
import pandas as pd

import tRep
import tRep.controller

__version__ = tRep.__version__

def main(**args):
    b6_loc = args.get('b6_loc')
    out_base = args.get('out_loc')
    database = args.get('database')

    # Load the hits
    Bdb = tRep.load_b6(b6_loc)

    # Load translations
    start = time.time()
    r2t = funcional_annotation(database, set(Bdb['target'].tolist()))
    end = time.time()
    print("{0:.1f} seconds to load translation database".format(end-start))

    # Make the translation table
    Bdb['functional_annotation'] = Bdb['target'].map(r2t)

    # Save
    Bdb.to_csv(os.path.join(out_base) + '_geneFunctionalAnnotation.tsv', \
        index=False, sep='\t')

def load_transtable(loc):
    db = pd.read_csv(loc, sep = '\t', names=['ID', 'annotation'], dtype=str)
    return db.set_index('ID')['annotation'].to_dict()

def funcional_annotation(tt_loc, hits):
    r2t = dict()
    assert type(hits) == type(set(hits))

    with gzip.open(tt_loc,'rt') as f:
        for line in f:
            if line.split('\t')[0] in hits:
                linewords = line.strip().split('\t')
                linewords[0] = linewords[1]
    return r2t


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple script to get functional taxonomy from a b6 file and a translation table")

    InpArgs = parser.add_argument_group('ARGUMENTS')
    InpArgs.add_argument('-b', '--b6_loc', help='location of b6+ file', \
        required=True)
    InpArgs.add_argument('-o', '--out_loc',  help=\
        'location of output file', required=True)
    InpArgs.add_argument('-d', '--database',  help=\
        'location of translation database', required=True)

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
