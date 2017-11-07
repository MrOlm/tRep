#!/usr/bin/env python

import pandas as pd

import tRep

def convert_b6_to_Tdb(args):
    b6_loc = args.get('b6_loc')
    out_loc = args.get('out_loc')

    #print("converting {0} to {1}".format(b6_loc, out_loc))

    # Load Bdb
    Bdb = tRep.load_b6(b6_loc)
    Bdb = Bdb[Bdb['taxID'].astype(str) != 'NA']

    # Add the taxonomy
    tax = tRep.gen_levels_db(list(Bdb['taxID'].unique()))
    Tdb = pd.merge(Bdb, tax, on='taxID', how='outer')

    # Save
    Tdb.to_csv(out_loc, index=False)
