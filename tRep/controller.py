#!/usr/bin/env python

import pandas as pd

import tRep

def convert_b6_to_Tdb(args, save=False):
    b6_loc = args.get('b6_loc')
    out_loc = args.get('out_loc')
    aa_loc = args.get('amino_acids', None)
    tax_type = args.get('tax_type', 'species')

    # Load Bdb
    Bdb = tRep.load_b6(b6_loc, tax_type=tax_type)
    Bdb = Bdb[~Bdb['taxID'].isna()]

    # Add the taxonomy
    tax = tRep.gen_levels_db(list(Bdb['taxID'].unique()))

    # Merge in
    Tdb = pd.merge(Bdb, tax, on='taxID', how='outer')
    assert len(Tdb) == len(Bdb)

    # Add back missing ones
    if aa_loc is not None:
        Adb = tRep.parse_prodigal_genes(aa_loc)

        # Figure out the overlap
        all_genes = set(Adb['gene'].tolist())
        hit_genes = set(Bdb['querry'].tolist())
        print("{0} of {1} genes have a {2} hit".format(len(hit_genes), len(all_genes), tax_type))
        assert len(hit_genes - all_genes) == 0

        # Add back some blanks
        db = pd.DataFrame({"querry":list(all_genes - hit_genes)})
        for level in tRep.get_levels():
            db[level] = 'unk'
        Tdb = pd.concat([Tdb, db]).reset_index(drop=True)

    # Save
    if save:
        Tdb.to_csv(out_loc, index=False)

    return Tdb

def gen_blank_levels():
    pass
