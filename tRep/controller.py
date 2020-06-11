#!/usr/bin/env python

import pandas as pd

import tRep

def extract_diamond_scaffold(id):
    '''
    From the gene name in diamond, extract the scaffold ID

    This is hard because diamond converts spaces in fasta headers to "|" symbols.
    '''
    # Try the way that works if | is not in the scaffold IDs
    scaff_gene = id.split('|')[0]

    # Make sure theres a gene available
    if (('_'  in scaff_gene) and (scaff_gene.split('_')[-1].isnumeric())):
        return '_'.join(scaff_gene.split('_')[:-1])

    # Otherwise try and parse as a non-diamond
    scaff_gene = id

    if (('_'  in scaff_gene) and (scaff_gene.split('_')[-1].isnumeric())):
        return '_'.join(scaff_gene.split('_')[:-1])

    # If you reach here, you dont know what to do
    print("I dont know how to parse the gene {0} into a scaffold; it will be ignored for scaffold and genome profiling".format(
        id))

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

    if len(Tdb) > 0:
        type = tRep.type_b6(b6_loc)
        if type == 'diamond':
            Tdb['scaffold'] = [extract_diamond_scaffold(x) for x in Tdb['querry']]
        elif type == 'b6+':
            Tdb['scaffold'] = ['_'.join(x.split('_')[:-1]) for x in Tdb['querry']]


    # Save
    if save:
        Tdb.to_csv(out_loc, index=False)

    return Tdb

def gen_blank_levels():
    pass
