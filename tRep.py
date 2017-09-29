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

from ete3 import NCBITaxa
ncbi = NCBITaxa()

import drep.d_cluster
import drep.d_bonus
import drep.d_filter

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
        tax = gen_taxonomy_string(Tdb[Tdb['score'] > 250], minPerc= int(\
            kwargs['percent']))
        print(tax)

    elif kwargs.get('method') == 'max':
        x = drep.d_bonus.gen_phylo_db(Tdb)
        taxID = x['tax_ID'][x['tax_confidence'] == x['tax_confidence'].max()].tolist()[0]
        tax = lineage_from_taxId(taxID)
        print(tax)

def gen_taxonomy_string(hits, minPerc= 50, testing=False):
    '''
    Determines the lowest taxonomic level with at least minPerc certainty

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
    for t in hits['taxID'].tolist():
        if t == 0:
            continue

        lin = ncbi.get_lineage(t)

        lin2name  = ncbi.get_taxid_translator(lin)
        name2rank = ncbi.get_rank(lin)

        #gen = False

        for i in lin:
            rank = name2rank[i]
            name = lin2name[i]
            if rank in countDic:
                #countDic[rank][name] = countDic[rank].get(name,0) + 1
                countDic[rank][i] = countDic[rank].get(i,0) + 1
            #if rank == 'phylum':
            #    gen = True

        #if gen == False:
        #    print(t)

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
