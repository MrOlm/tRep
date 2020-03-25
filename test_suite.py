#!/usr/bin/env python
'''
Run test of tRep
'''

import os
import glob
import shutil

import pandas as pd
from subprocess import call

import tRep

def load_b6_loc():
    return os.path.join(str(os.getcwd()), \
        'testFiles/N1_003_000G1_scaffold_min1000.fa.genes.faa-vs-uniprot.b6+')

def load_testdir():
    return os.path.join(str(os.getcwd()), \
        'testFiles/')

def load_ggkbase_loc():
    return os.path.join(str(os.getcwd()), \
        'testFiles/N1_003_000G1.contig-taxonomy.tsv')

def load_ggkbase_org_table():
    return os.path.join(str(os.getcwd()), \
        'testFiles/N1_003_000G1.organism_info.tsv')

def load_sample_Tdb():
    return os.path.join(str(os.getcwd()), \
        'testFiles/Tdb.csv')

def load_sample_s2b_loc():
    return os.path.join(str(os.getcwd()), \
        'testFiles/N1_003.delta.stb')

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'testFiles/test_backend/testdir/')
    return loc

def get_script_loc(script):
    if script == 'tax_collector.py':
        return os.path.join(str(os.getcwd()), \
            'bin/tax_collector.py')
    elif script == 'make_Tdb.py':
        return os.path.join(str(os.getcwd()), \
            'bin/make_Tdb.py')
    elif script == 'functional_tax.py':
        return os.path.join(str(os.getcwd()), \
            'bin/functional_tax.py')

class testUniProt():
    def setUp(self):
        self.b6_loc = load_b6_loc()
        self.ggkbase_s2t_loc = load_ggkbase_loc()

    def tearDown(self):
        pass

    def run(self):
        self.setUp()
        self.main_test_1()
        self.main_test_2()
        self.main_test_3()
        self.main_test_4()
        self.tearDown()

    def main_test_1(self):
        '''
        test loading b6 file and running on large amounts of data
        '''
        Bdb = tRep.load_b6(self.b6_loc)
        Bdb = Bdb[~Bdb['taxID'].isna()]
        tax = tRep.gen_taxonomy_string(Bdb['taxID'].tolist())
        assert tax == '1239|Bacteria|Firmicutes|unk|unk|unk|unk|unk'

    def main_test_2(self):
        '''
        test gen_levels_db (making Tdb)
        '''
        # Load Bdb
        Bdb = tRep.load_b6(self.b6_loc)
        Bdb = Bdb[~Bdb['taxID'].isna()]

        # Get Tdb
        tax = tRep.gen_levels_db(list(Bdb['taxID'].unique()))
        Tdb = pd.merge(Bdb, tax, on='taxID', how='outer')

        assert len(Tdb) == len(Bdb), [len(Tdb), len(Bdb)]

    def main_test_3(self):
        '''
        Make scaffold level taxonomy calls; compare to ggkbase taxonomy calls
        '''
        Tdb = pd.read_csv(load_sample_Tdb())
        Gdb = pd.read_table(self.ggkbase_s2t_loc)

        # Generate scaffold2taxonomy
        Sdb = tRep.gen_taxonomy_table(Tdb, on='scaffold')

        # compare genus-level calls = make sure at least half are the same
        for level in ['genus']:#, 'species']:
            s_col = level + '_winner'
            g_col = level[0].upper() + level[1:] + ' winner'

            s2t = Sdb.set_index('scaffold')[s_col].to_dict()
            Gdb[s_col] = Gdb['Contig name'].map(s2t)
            tdb = Gdb[Gdb[s_col] != Gdb[g_col]]
            assert len(tdb)/len(Gdb) < .5

    def main_test_4(self):
        '''
        Test generating taxonomy based on scaffold 2 bin file
        '''
        Tdb = pd.read_csv(load_sample_Tdb())
        stb = tRep.load_stb(load_sample_s2b_loc())

        # Add the bin to Tdb
        Tdb = tRep.add_bin_to_tdb(Tdb, load_sample_s2b_loc())
        Tdb['bin'] = ['N1_003_000G1_UNK' if x == 'unk' else x.replace('.', '_') \
            for x in Tdb['bin']]

        # Load ggkbase org table
        Gdb = pd.read_table(load_ggkbase_org_table())

        # Make sure the ggkbase org table and my table have the same bins
        assert len(Gdb['name'].unique()) == len(Tdb['bin'].unique())
        for org in Gdb['name'].unique():
            assert org in Tdb['bin'].tolist(), org

        # Make Sdb
        Sdb = tRep.gen_taxonomy_table(Tdb, on='bin')

        # Make sure Sdb and ggkbase get the same taxonomy
        Gdb['new_tax'] = Gdb['name'].map(Sdb.set_index('bin')['taxonomy'].to_dict())
        db = Gdb[['new_tax', 'name', 'taxonomy']]
        for i, row in db.iterrows():
            if (row['name'] in ['N1_003_000G1_UNK', 'dasN1_003_000G1_concoct_26_fa']):
                continue

            my_tax = row['new_tax'].split('|')[-1]
            gg_tax = row['taxonomy'].split(',')[0]

            assert my_tax == gg_tax

class testMakeTdb():
    def setUp(self):
        self.b6_loc = load_b6_loc()
        self.diamond_loc = load_testdir() + 'N4_005_008G1_Pseudomonas_aeruginosa_66_425.proteins.translated.diamondOut'
        self.test_dir = load_random_test_dir()
        self.script_loc = get_script_loc('make_Tdb.py')

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.main_test_1()
        self.tearDown()

        self.setUp()
        self.main_test_2()
        self.tearDown()

    def main_test_1(self):
        '''
        Make sure MakeTdb makes all the files when run by default
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.b6_loc, '-o', out_base]
        call(cmd)

        files = glob.glob(out_base + '*')
        assert len(files) == 1
        for f in files:
            assert os.path.getsize(f) > 0

    def main_test_2(self):
        '''
        Make sure MakeTdb makes all the files when run with diamondOut
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base]
        call(cmd)

        files = glob.glob(out_base + '*')
        print(files)
        assert len(files) == 1
        for f in files:
            assert os.path.getsize(f) > 0

class testFunctional_tax():
    def setUp(self):
        self.b6_loc = load_b6_loc()
        self.test_dir = load_random_test_dir()
        self.diamond_loc = load_testdir() + 'N4_005_008G1_Pseudomonas_aeruginosa_66_425.proteins.translated.diamondOut'
        self.aa_loc = load_testdir() + 'N4_005_008G1_Pseudomonas_aeruginosa_66_425.proteins.faa'
        self.script_loc = get_script_loc('functional_tax.py')
        self.small_tans_table = load_testdir() + 'uniref100.subset.ttable.gz'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test_1()
        self.tearDown()

    def test_1(self):
        out_base = os.path.join(self.test_dir, 'test_out_base')
        cmd = [self.script_loc, '-b', self.b6_loc, '-o', out_base, '-d', self.small_tans_table]
        call(cmd)

        files = glob.glob(out_base + '*')
        assert len(files) == 1, files
        print(files)

class test_tax_collector():
    def setUp(self):
        self.b6_loc = load_b6_loc()
        self.test_dir = load_random_test_dir()
        self.diamond_loc = load_testdir() + 'N4_005_008G1_Pseudomonas_aeruginosa_66_425.proteins.translated.diamondOut'
        self.diamond_loc_nospecies = load_testdir() + 'no_species.diamondOut'
        self.aa_loc = load_testdir() + 'N4_005_008G1_Pseudomonas_aeruginosa_66_425.proteins.faa'
        self.script_loc = get_script_loc('tax_collector.py')
        self.ggkbase_s2t_loc = load_ggkbase_loc()

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.main_test_1()
        self.tearDown()

        self.setUp()
        self.main_test_2()
        self.tearDown()

        self.setUp()
        self.main_test_3()
        self.tearDown()

        self.setUp()
        self.main_test_4()
        self.tearDown()

        self.setUp()
        self.main_test_5()
        self.tearDown()

        # THIS TAKES A WHILE AND REQUIRES SUPERVISION; NO NEED TO RUN NORMALLY
        # self.setUp()
        # self.update_test()
        # self.tearDown()

    def main_test_1(self):
        '''
        Make sure tax collector makes all the files when run by default
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.b6_loc, '-o', out_base, '-stb', 'ALL']
        call(cmd)

        files = glob.glob(out_base + '*')
        assert len(files) == 3
        for f in files:
            if 'fullGenomeTaxonomy' in f:
                db = pd.read_csv(f)
                assert len(db) == 1
            assert os.path.getsize(f) > 0

        # Do it with an .stb file
        cmd = [self.script_loc, '-b', self.b6_loc, '-o', out_base, '-stb', load_sample_s2b_loc(), '--SkipScaffolds']
        call(cmd)

        files = glob.glob(out_base + '*')
        for f in files:
            if 'fullGenomeTaxonomy' in f:
                db = pd.read_csv(f, sep='\t')
                assert len(db) > 1
            assert os.path.getsize(f) > 0

    def main_test_2(self):
        '''
        Make sure tax collector makes dimaond when run by default
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base]
        call(cmd)

        files = glob.glob(out_base + '*')
        assert len(files) == 2
        for f in files:
            assert os.path.getsize(f) > 0

    def main_test_3(self):
        '''
        Make sure tax collector works with aa
        '''
        # Try with group-level taxonomy
        out_base = os.path.join(self.test_dir, 'test_out_base_g')
        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base, '-a', self.aa_loc, '--tax_type', 'group', '-stb', 'ALL']
        call(cmd)

        # Load gene taxonomy
        file = glob.glob(out_base + '*fullGeneTaxonomy.tsv')[0]
        gdb = pd.read_csv(file, sep='\t')
        assert len(gdb) == 6282, len(gdb)

        # Load scaffold taxonomy
        file = glob.glob(out_base + '*fullScaffoldTaxonomy.tsv')[0]
        sdb = pd.read_csv(file, sep='\t')
        assert len(sdb) == 59

        # Load genome taxonomy
        file = glob.glob(out_base + '*fullGenomeTaxonomy.tsv')[0]
        genomedb = pd.read_csv(file, sep='\t')
        assert len(genomedb) == 1
        print("Group - {0}".format(genomedb['full_taxonomy'].tolist()[0]))

        # Try with species-level taxonomy
        out_base = os.path.join(self.test_dir, 'test_out_base_s')
        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base, '-a', self.aa_loc, '--tax_type', 'species', '-stb', 'ALL']
        call(cmd)

        # Load gene taxonomy
        file = glob.glob(out_base + '*fullGeneTaxonomy.tsv')[0]
        gdbS = pd.read_csv(file, sep='\t')
        assert len(gdbS) == 6282, len(gdb)

        # Load scaffold taxonomy
        file = glob.glob(out_base + '*fullScaffoldTaxonomy.tsv')[0]
        sdbS = pd.read_csv(file, sep='\t')
        assert len(sdbS) == 59

        # Load genome taxonomy
        file = glob.glob(out_base + '*fullGenomeTaxonomy.tsv')[0]
        genomedbS = pd.read_csv(file, sep='\t')
        assert len(genomedbS) == 1
        print("Species - {0}".format(genomedbS['full_taxonomy'].tolist()[0]))

        # Try without adding aa
        out_base = os.path.join(self.test_dir, 'test_out_base_N')
        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base, '--tax_type', 'species', '-stb', 'ALL']
        call(cmd)

        # Load gene taxonomy
        file = glob.glob(out_base + '*fullGeneTaxonomy.tsv')[0]
        gdbN = pd.read_csv(file, sep='\t')
        assert len(gdbN) != 6282, len(gdb)

        # Load scaffold taxonomy
        file = glob.glob(out_base + '*fullScaffoldTaxonomy.tsv')[0]
        sdbN = pd.read_csv(file, sep='\t')
        assert len(sdbN) != 59, len(sdbN)

        # Load genome taxonomy
        file = glob.glob(out_base + '*fullGenomeTaxonomy.tsv')[0]
        genomedbN = pd.read_csv(file, sep='\t')
        assert len(genomedbN) == 1

    def main_test_4(self):
        '''
        Make sure tax collector works when the taxID hit has no species call
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.diamond_loc_nospecies, '-o', out_base]
        call(cmd)

        # Load the gene table
        Gdb = pd.read_csv(glob.glob(out_base + '*')[0], sep='\t')
        assert Gdb['superkingdom_winner'].tolist()[0] == 'Bacteria'

        files = glob.glob(out_base + '*')
        assert len(files) == 2
        for f in files:
            assert os.path.getsize(f) > 0

    def main_test_5(self):
        '''
        Compare group and species taxonomy
        '''
        # Try with group-level taxonomy
        out_base = os.path.join(self.test_dir, 'test_out_base_g')
        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base, '-a', self.aa_loc, '--tax_type', 'group', '-stb', 'ALL']
        call(cmd)

        # Load gene taxonomy
        file = glob.glob(out_base + '*fullGeneTaxonomy.tsv')[0]
        ggdb = pd.read_csv(file, sep='\t')
        assert len(ggdb) == 6282, len(ggdb)

        # Try with species-level taxonomy
        out_base = os.path.join(self.test_dir, 'test_out_base_s')
        cmd = [self.script_loc, '-b', self.diamond_loc, '-o', out_base, '-a', self.aa_loc, '--tax_type', 'species', '-stb', 'ALL']
        call(cmd)

        # Load gene taxonomy
        file = glob.glob(out_base + '*fullGeneTaxonomy.tsv')[0]
        gsdb = pd.read_csv(file, sep='\t')
        assert len(gsdb) == 6282, len(gsdb)

        # Compare
        # Find cases where species is missing
        missing_species = set(gsdb[gsdb['superkingdom'] == 'unk']['querry'].tolist())

        # Get the group hits
        db = ggdb[(ggdb['querry'].isin(missing_species)) & (ggdb['superkingdom'] != 'unk')]
        db['problem_species'] = [x.split('_')[1] if len(x.split('_')) > 2 else 'WEIRD' for x in db['target']]
        print(db['problem_species'].value_counts())

    def update_test(self):
        '''
        Make sure tax collector actually updates
        '''
        out_base = os.path.join(self.test_dir, 'test_out_base')

        cmd = [self.script_loc, '-b', self.diamond_loc_nospecies, '-o', out_base, '--update']
        call(cmd)

        # Load the gene table
        Gdb = pd.read_csv(glob.glob(out_base + '*')[0], sep='\t')
        assert Gdb['superkingdom_winner'].tolist()[0] == 'Bacteria'

        files = glob.glob(out_base + '*')
        assert len(files) == 2
        for f in files:
            assert os.path.getsize(f) > 0

if __name__ == '__main__':
    testUniProt().run()
    test_tax_collector().run()
    testMakeTdb().run()
    testFunctional_tax().run()
    print('everything is working swimmingly!')
