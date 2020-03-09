#!/usr/bin/env python

"""
tRep
"""
import os
__author__ = "Matt Olm"
__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                'VERSION')).read().strip()
__license__ = "MIT"

import glob
import time
import shutil
import argparse
import tempfile
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

from Bio import SeqIO
from ete3 import NCBITaxa
from collections import defaultdict
ncbi = NCBITaxa()

import drep.d_bonus
import drep.d_filter
import drep.d_cluster

'''
THIS SECTION IS BASED ON USEARCH
'''
def type_b6(location):
    '''
    Figure out the type of b6+ file this actually is

    All based on the number of underscores is the second column
    '''
    types = []
    count = 0
    with open(location, 'r') as o:
        for line in o.readlines():
            line = line.strip()
            underscores = line.split()[1].count('_')
            types.append(underscores)
            count += 1
            if count >= 5:
                break

    if set(types) == set([1]):
        return 'b6+'
    elif set(types) == set([2]):
        return 'diamond'
    else:
        print("I cant tell what kind of b6+ file you have! Quitting")
        raise Exception()

def load_b6(location, tax_type='species'):
    '''
    return the b6 file as a pandas DataFrame
    '''
    type = type_b6(location)

    if type == 'b6+':
        '''
        parse the b6+ format, like from Brian
        '''
        header = ['querry', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',\
            'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',\
            'extra']
        Bdb = pd.read_csv(location, names=header, sep='\t')
        Bdb['annotation'], Bdb['taxID'], Bdb['taxString'] = zip(*Bdb['extra'].map(parse_b6))
        Bdb['scaffold'] = ['_'.join(x.split('_')[:-1]) for x in Bdb['querry']]

    elif type == 'diamond':
        '''
        parse like this is from the diamond out
        '''
        header = ['querry', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',\
            'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',\
            'extra']
        Bdb = pd.read_csv(location, names=header, sep='\t')
        Bdb['taxID'] = Bdb['target'].apply(parse_diamond, tax_type=tax_type)
        Bdb['scaffold'] = ['_'.join(x.split('|')[0].split('_')[:-1]) for x in Bdb['querry']]

    else:
        print("I dont know how to parse type {0}".format(type))
        raise Exception()

    return Bdb

def parse_prodigal_genes(gene_fasta):
    '''
    Parse the prodigal .fna file

    Return a datatable with gene info and a dictionary of gene -> sequence
    '''
    table = defaultdict(list)
    print(gene_fasta)
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)

        table['gene'].append(gene)
        table['scaffold'].append("_".join(gene.split("_")[:-1]))
        table['direction'].append(record.description.split("#")[3].strip())
        table['partial'].append('partial=00' not in record.description)

        # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE
        table['start'].append(int(record.description.split("#")[1].strip())-1)
        table['end'].append(int(record.description.split("#")[2].strip())-1)

    Gdb = pd.DataFrame(table)

    return Gdb

def parse_diamond(line, tax_type='species'):
    if tax_type == 'species':
        loc = 1
    elif tax_type == 'group':
        loc = 2
    try:
        taxID = float(line.split('_')[loc])
    except:
        taxID = np.nan
    return taxID

def parse_b6(line):
    '''
    return a parsed list from a b6+ annotation string:

    [annotation, \ # the functional annotation
    taxID, \
    taxString
    ]
    taxOrder = ['species', 'phyla', 'class', 'order', 'family', 'genus']
    '''
    words = [x.strip() for x in line.split(';')]
    annotation = words[0]
    try:
        taxID = float([x for x in words if x.startswith('TaxID')][0].split()[0].replace(\
                'TaxID=',''))
    except:
        taxID = np.nan

    try:
        taxString = line.split('"')[1]
        taxList = [x.strip() for x in taxString.split(';')]
    except:
        taxString = ';'.join(['NA'] * 7)

    return annotation, taxID, taxString

def gen_levels_db(hits):
    '''
    from a list of taxIDs, return a DataFrame deliniating their taxonomies
    '''
    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # start the spinner
    # spinner = Halo(text='Parsing taxIDs', spinner='dots')
    # spinner.start()
    print("parsing taxIDs...")

    # set up the database
    table = {'taxID':[]}
    for level in Levels:
        table[level] = []

    # fill in nested dictionary
    for t in hits:
        if t == 0:
            continue

        try:
            lin = ncbi.get_lineage(t)
        except:
            continue

        lin2name  = ncbi.get_taxid_translator(lin)
        name2rank = ncbi.get_rank(lin)
        rank2name = {v: k for k, v in name2rank.items()}

        for level in Levels:
            if level in rank2name:
                table[level].append(lin2name[rank2name[level]])
            else:
                table[level].append('unk')
        table['taxID'].append(t)

    # spinner.stop()

    return pd.DataFrame(table)

def get_levels():
    return ['superkingdom','phylum','class','order','family','genus','species']

def gen_taxonomy_table(Idb, on='scaffold', minPerc=50):
    '''
    From a dataframe with all of the levels present, calculate percentages
    '''
    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # # start the spinner
    # spinner = Halo(text='Generating taxonomy table', spinner='dots')
    # spinner.start()
    print("Generating taxonomy table...")

    # set up Sdb
    table = {}
    for level in Levels:
        table[level + '_winner'] = []
        table[level + '_percent'] = []
    table[on] = []

    # make Sdb
    for thing, db in Idb.groupby(on):
        table[on].append(thing)
        for level in Levels:
            vcounts = dict(db[level].value_counts())
            try:
                winner = max(vcounts, key=vcounts.get)
            except:
                table[level + '_winner'].append('unk')
                table[level + '_percent'].append(100)
                continue

            total = sum(vcounts.values())
            count = vcounts[winner]
            table[level + '_winner'].append(winner)
            table[level + '_percent'].append(((count/total) *100))
    Sdb = pd.DataFrame(table)

    # add taxonomy
    Sdb['minPerc'] = int(minPerc)
    Sdb['full_taxonomy'] = Sdb.apply(calculate_full_taxonmy, axis=1)
    Sdb['taxonomy'] = Sdb['full_taxonomy'].map(get_simple_tax)
    del Sdb['minPerc']

    # spinner.stop()
    return Sdb

def get_simple_tax(full_taxonomy):
    '''
    from the full taxonomy, get the simple taxonomy (just the highest name)
    '''
    taxes = full_taxonomy.split('|')
    for t in reversed(taxes):
        if t != 'unk':
            return t
    return 'unk'

def calculate_full_taxonmy(row):
    '''
    calculate consensus taxonomy based on a row of Sdb
    '''
    minPerc = row['minPerc']
    Levels = ['superkingdom','phylum','class','order','family','genus','species']
    taxString = ''
    skip = False
    for l in Levels:
        if (int(row[l + '_percent']) < minPerc) | skip:
            taxString += 'unk|'
            skip = True
        else:
            taxString += row[l + '_winner'] + '|'
    return taxString[:-1]

def add_bin_to_tdb(tdb, s2b):
    '''
    Add the 'bin' column to tdb based off of an stb file
    '''
    if type(s2b) == dict:
        pass
    else:
        s2b = load_stb(s2b)

    # Add the bin to Tdb
    tdb['bin'] = tdb['scaffold'].map(s2b)
    tdb['bin'] = tdb['bin'].fillna('unk')

    # Make sure at least some scaffolds map
    assert len(tdb[tdb['bin'] == 'unk']) != len(tdb), \
        "No scaffolds in the stb map to a bin"
    assert len(tdb[tdb['bin'] != tdb['bin']]) == 0

    return tdb

def load_stb(file):
    stb = {}
    with open(file,'r') as ins:
        for line in ins:
            words = line.strip().split('\t')
            if words[0].startswith('scaffold_name'): continue
            scaffold = words[0]
            b = words[1]
            stb[scaffold] = b
    return stb

###############################################################################
'''
THIS SECTION IS BASED ON CENTRIFUGE
'''

def main(**kwargs):
    """ Main entry point of the app """
    cent_hit = kwargs.get('cent_hit')
    cent_report = kwargs.get('cent_report')
    fasta = kwargs.pop('fasta')

    if (cent_hit != None) and (cent_report != None):
        print('determining taxonomy from raw centrifuge')
        from_raw_centrifuge(cent_hit, cent_report, **kwargs)

    elif fasta != None:
        print('determining taxonomy from .fasta file')
        from_fasta(fasta, **kwargs)

def from_fasta(fasta, **kwargs):
    dirpath = tempfile.mkdtemp() + '/'

    bdb = pd.DataFrame({'location':[fasta], 'genome':[os.path.basename(fasta)]})
    print('running prodigal')
    drep.d_filter.run_prodigal(bdb, dirpath)
    print('running centrifuge')
    drep.d_bonus.run_centrifuge(bdb, dirpath, dirpath, cent_index=kwargs['cent_index'])

    hits, report = glob.glob(dirpath + os.path.basename(fasta) + '*.tsv')
    from_raw_centrifuge(hits, report, **kwargs)

def from_raw_centrifuge(hits, report, **kwargs):
    '''
    Return the taxonomy based on raw centrifuge output
    '''
    Tdb = drep.d_bonus.parse_raw_centrifuge(hits, report)
    Tdb = filter_hits(Tdb, **kwargs)
    from_tdb(Tdb, **kwargs)

def filter_hits(Tdb, **kwargs):
    min_score = int(kwargs.get('min_score', 250))
    min_diff = int(kwargs.get('min_diff', 0))

    Tdb['diff'] = [x - y for x, y in zip(Tdb['score'], Tdb['2ndBestScore'])]

    Tdb = Tdb[Tdb['score'] >= min_score]
    Tdb = Tdb[Tdb['diff'] >= min_diff]

    return Tdb

def from_tdb(Tdb, **kwargs):
    if kwargs.get('full_dump') != None:
        tdb = gen_full_tdb(Tdb)
        tdb.to_csv(kwargs.get('full_dump'), index=False)

    elif kwargs.get('method') == 'percent':
        tax = gen_taxonomy_string(Tdb['hit'][Tdb['score'] > 250], minPerc= int(\
            kwargs['percent']).to_list())
        print(tax)

    elif kwargs.get('method') == 'max':
        x = drep.d_bonus.gen_phylo_db(Tdb)
        taxID = x['tax_ID'][x['tax_confidence'] == x['tax_confidence'].max()].tolist()[0]
        tax = lineage_from_taxId(taxID)
        print(tax)

def gen_taxonomy_string(hits, minPerc= 50, testing=False):
    '''
    Determines the lowest taxonomic level with at least minPerc certainty

    hits = list of taxIDs

    For every hit:
        reconstruct the lineage (kingdom, phylum, class, ect.)
        add a count to every rank in the lineage

    For every rank:
        see if the number of hits matching one taxa at that rank is above the minPerc
        the denominator for this equation is the number of hits that have a phyla rank

    * Note: this is complicated because some lower ranks don't have higher ranks
        For example, species [Eubacterium] rectale (taxID 39491) has no genus
        Also, species [artifical construct] (taxID 32630) has no anything but species

    '''
    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # generate nested dictionary for levels
    countDic = {}
    for level in Levels:
        countDic[level] = {}

    # fill in nested dictionary
    for t in hits:
        if t == 0:
            continue

        try:
            lin = ncbi.get_lineage(t)
        except:
            continue

        lin2name  = ncbi.get_taxid_translator(lin)
        name2rank = ncbi.get_rank(lin)

        for i in lin:
            rank = name2rank[i]
            name = lin2name[i]
            if rank in countDic:
                countDic[rank][i] = countDic[rank].get(i,0) + 1

    # make the taxonomy string
    WinningName = None

    # figure out the total number of hits
    total = sum(countDic['phylum'].values())

    # find the lowest rank winner
    winner = False

    Levels.reverse()
    for level in Levels:
        if winner != False:
            break
        dic = countDic[level]
        for name in sorted(dic, key=dic.get, reverse= True):
            count = dic[name]
            if ((count/total) *100) > minPerc:
                winner = name
            break

    # return the whole taxonomy string
    name = []

    lin = ncbi.get_lineage(winner)

    lin2name  = ncbi.get_taxid_translator(lin)
    name2rank = ncbi.get_rank(lin)
    rank2name = {v: k for k, v in name2rank.items()}

    Levels.reverse()
    for level in Levels:
        if level in rank2name:
            name.append(lin2name[rank2name[level]])
        else:
            name.append('unk')

    if testing:
        return '|'.join([str(winner)] + name), countDic
    else:
        return '|'.join([str(winner)] + name)

def gen_full_tdb(hits):
    '''
    Report all hits to all taxonomic levels
    '''

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # generate nested dictionary for levels
    countDic = {}
    for level in Levels:
        countDic[level] = {}

    # fill in nested dictionary
    for t in hits['taxID'].tolist():
        if t == 0:
            continue

        # This try / except thing is trying to catch sporatic errors of:
        # sqlite3.OperationalError: disk I/O error
        try:
            lin = ncbi.get_lineage(t)
            lin2name  = ncbi.get_taxid_translator(lin)
            name2rank = ncbi.get_rank(lin)
        except:
            time.sleep(1)
            lin = ncbi.get_lineage(t)
            lin2name  = ncbi.get_taxid_translator(lin)
            name2rank = ncbi.get_rank(lin)

        printed=False
        for j, i in enumerate(lin[::-1]):
            rank = name2rank[i]
            name = lin2name[i]

            if rank in countDic:
                if printed == False:
                    #print(rank)
                    printed=True

                countDic[rank][name] = countDic[rank].get(name,0) + 1

    # make the table
    total = sum(countDic['phylum'].values())
    table = {'tax_confidence':[], 'tax_level':[], 'taxonomy':[]}

    for level in Levels:
        dic = countDic[level]
        for name in sorted(dic, key=dic.get, reverse= True):
            count = dic[name]

            table['tax_confidence'].append(((count/total) *100))
            table['tax_level'].append(level)
            table['taxonomy'].append(name)

    tdb = pd.DataFrame(table)
    return tdb

def lineage_from_taxId(t):
    Levels = ['superkingdom','phylum','class','order','family','genus','species']
    name = []

    lin = ncbi.get_lineage(t)

    lin2name  = ncbi.get_taxid_translator(lin)
    name2rank = ncbi.get_rank(lin)
    rank2name = {v: k for k, v in name2rank.items()}

    for level in Levels:
        if level in rank2name:
            name.append(lin2name[rank2name[level]])
        else:
            name.append('unk')

    return '|'.join([str(int(t))] + name)
