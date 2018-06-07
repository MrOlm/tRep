#!/usr/bin/env python

import pandas as pd

import tRep

def convert_b6_to_Tdb(args, save=False):
    b6_loc = args.get('b6_loc')
    out_loc = args.get('out_loc')

    # Load Bdb
    Bdb = tRep.load_b6(b6_loc)
    Bdb = Bdb[~Bdb['taxID'].isna()]

    # Add the taxonomy
    tax = tRep.gen_levels_db(list(Bdb['taxID'].unique()))

    # Merge in
    Tdb = pd.merge(Bdb, tax, on='taxID', how='outer')

    assert len(Tdb) == len(Bdb)

    # Save
    if save:
        Tdb.to_csv(out_loc, index=False)

    return Tdb
