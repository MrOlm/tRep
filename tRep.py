#!/usr/bin/env python

"""
tRep
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import tempfile
import shutil
import glob
import pandas as pd
import os
import time
from halo import Halo

import warnings
warnings.filterwarnings("ignore")

from ete3 import NCBITaxa
ncbi = NCBITaxa()

import drep.d_cluster
import drep.d_bonus
import drep.d_filter

'''
THIS SECTION IS BASED ON USEARCH
'''
def load_b6(location, type='b6+'):
    '''
    return the b6 file as a pandas DataFrame
    '''

    Bdb = pd.read_table(location, header=None)
    if type == 'b6+':
        '''
        parse the b6+ format, like from Brian
        '''
        header = ['querry', 'target', 'percentID', 'alignment_length', 'mm', 'gaps',\
            'querry_start', 'querry_end', 'target_start', 'target_end', 'e-value', 'bit_score',\
            'extra']
        Bdb.columns = header
        Bdb['annotation'], Bdb['taxID'], Bdb['taxString'] = zip(*Bdb['extra'].map(parse_b6))
        Bdb['scaffold'] = ['_'.join(x.split('_')[:-1]) for x in Bdb['querry']]

    return Bdb

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
        taxID = int([x for x in words if x.startswith('TaxID')][0].split()[0].replace(\
                'TaxID=',''))
    except:
        taxID = 'NA'

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
    spinner = Halo(text='Parsing taxIDs', spinner='dots')
    spinner.start()

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

    spinner.stop()

    return pd.DataFrame(table)

def gen_taxonomy_table(Idb, on='scaffold', minPerc=50):
    '''
    From a dataframe with all of the levels present, calculate percentages
    '''
    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # start the spinner
    spinner = Halo(text='Generating taxonomy table', spinner='dots')
    spinner.start()

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

    spinner.stop()
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

    args = parser.parse_args()
    main(**vars(args))
