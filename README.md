# tRep
Some scripts to help you get the taxonomy of a genome

### Overview
This program is a collection of tools for determining the taxonomy of a genome using a number of different tools. Its pretty unorgaized at the moment, but might develop into something at some point.

### Installation

```
pip install taxRep
```

### Dependencies

dRep (https://github.com/MrOlm/drep)

ETE3 (https://github.com/etetoolkit/ete)

## Data preparation

At the heart of much of tRep's functionality is the so-called b6+ file. This file is the result of performing an amino acid search of the sequences you're interested in against the UniRef100 database.There are two ways of making this file, but neither of them are easy at this point.

### Annotation using usearch
If you're in the Banfield lab, you make these during the standard processing pipeline using usearch, with a command that looks like:

```
/opt/bin/bio/usearch64 -ublast ${file} -db /data1/blastdb/uniref/uniprot_cp.fasta.udb -maxhits 1 -evalue 0.0001 -threads 48 -blast6out ${file}.b6.fix

ruby /home/bcthomas/work/pipeline_dev/fix_usearch_output.rb ${file}.b6.fix > $file.b6

annolookup.py ${file}.b6 uniprot > $file.b6+
```

However, this requires 64bit usearch because the database is so large (which costs lots of money) and the script annolookup.py to translate from UniRef100 ID to NCBI TaxID (but the script relies on a special server, Hydra, which is not public). In the end you get a b6+ file that looks like:

```
N1_003_000G1_scaffold_0_38	A0A0F4WVY9_9CLOT	100.0	55	0	0	1	55	1	55	2.9e-25	121.3	Carbonic anhydrase, family 3 {ECO:0000313|EMBL:KJZ97909.1};  TaxID=1523156 species="Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;" source="Clostridium sp. IBUN13A.;"
N1_003_000G1_scaffold_0_26	A0A0F4VQV2_9CLOT	100.0	83	0	0	1	83	1	83	3.5e-38	164.9	Uncharacterized protein {ECO:0000313|EMBL:KJZ83027.1};  TaxID=1523154 species="Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;" source="Clostridium sp. IBUN125C.;"
N1_003_000G1_scaffold_0_12	A0A0F4VP19_9CLOT	100.0	109	0	0	1	109	1	109	7e-55	220.7	Uncharacterized protein {ECO:0000313|EMBL:KJZ83035.1};  TaxID=1523154 species="Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;" source="Clostridium sp. IBUN125C.;"
N1_003_000G1_scaffold_0_20	C4IHT5_CLOBU	100.0	100	0	0	1	100	1	100	1.2e-48	199.9	Uncharacterized protein {ECO:0000313|EMBL:EEP53634.1};  TaxID=632245 species="Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;" source="Clostridium butyricum E4 str. BoNT E BL5262.;"
N1_003_000G1_scaffold_0_10	R0ASX2_CLOBU	100.0	62	0	0	1	62	1	62	3.3e-25	121.3	Uncharacterized protein {ECO:0000313|EMBL:ENZ35859.1};  TaxID=997898 species="Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;" source="Clostridium butyricum 60E.3.;"
```

### Annotation using diamond

A free alternative to usearch is the program Diamond. A way of avoiding a special server to translate between UniRef100 and NCBI TaxIDs is to simply put the TaxID in the diamond database (an idea I borrowed from the program Tango (https://github.com/NBISweden/tango)). This means you can then download this special diamond database, search against it, and immediately have a b6+ file.

If you're in the Sonnenburg lab, this special database is located at `/LAB_DATA/DATABASES/UniRef100/uniref100.translated.diamond.dmnd`, and has a companion table of annotations at `/LAB_DATA/DATABASES/UniRef100/uniref100.ttable`. The entries in this database look like:

`>Q6GZX2_2219561_10492`

The first part is the UniRef100 ID; the second part is the taxID of the representative organism of the cluster (based on the file `ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz`); the third part is the taxID of the UniRef100 cluster as a whole (based on the original fasta header description).

If you're not in the Sonnenburg lab and are interested in using this tool, drop me a line. I'm happy to share this database, but because it is ~80Gb, I can't think of an easy way to share it publicly. You can also recreate this yourself if you'd like using the format described above.

You can then search against this database using a command like:

```
diamond blastp -d  /LAB_DATA/DATABASES/UniRef100/uniref100.translated.diamond.dmnd -p 6 -q ${file} -o ${file}.diamondOut -f 6 -e 0.0001 -k 1
```

Which will make a b6+ file that looks like:

```
N4_005_008G1_scaffold_0_1|N4_005_008G1_Pseudomonas_aeruginosa_66_425|N4_005_008G1       A0A0F7QWA2_294_286      100.0   159     0       0       1       159     1       159     2.8e-83      317.8
N4_005_008G1_scaffold_0_2|N4_005_008G1_Pseudomonas_aeruginosa_66_425|N4_005_008G1       A0A0H2ZJ63_208963_287   100.0   638     0       0       1       638     1       638     0.0e+00      1226.5
N4_005_008G1_scaffold_0_3|N4_005_008G1_Pseudomonas_aeruginosa_66_425|N4_005_008G1       A0A263PVC1_1350465_286  100.0   192     0       0       1       192     1       192     5.7e-99      370.2
N4_005_008G1_scaffold_0_4|N4_005_008G1_Pseudomonas_aeruginosa_66_425|N4_005_008G1       A0A077K0S4_1350465_286  100.0   252     0       0       1       252     1       252     5.8e-128     466.8
N4_005_008G1_scaffold_0_5|N4_005_008G1_Pseudomonas_aeruginosa_66_425|N4_005_008G1       A0A072ZJM5_294_286      100.0   209     0       0       1       209     1       209     8.0e-107     396.4
```

Where the first column is the query sequence and the second column is the modified header with taxID information.

**Yes, I know both of these files look different, but tRep is able to parse them both just fine if you use them when it asks for a b6+ file**

## make_Tdb.py

This program makes a Tdb.csv file from a b6+ file resulting from an assembly. A Tdb.csv file contains lots of taxonomy information that can then be used by other aspects of the program to actually get the taxonomy that you want. A Tdb.csv file looks like this:

```
querry,target,percentID,alignment_length,mm,gaps,querry_start,querry_end,target_start,target_end,e-value,bit_score,extra,annotation,taxID,taxString,scaffold,class,family,genus,order,phylum,species,superkingdom
N1_003_000G1_scaffold_0_38,A0A0F4WVY9_9CLOT,100.0,55,0,0,1,55,1,55,2.9e-25,121.3,"Carbonic anhydrase, family 3 {ECO:0000313|EMBL:KJZ97909.1};  TaxID=1523156 species=""Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;"" source=""Clostridium sp. IBUN13A.;"" ","Carbonic anhydrase, family 3 {ECO:0000313|EMBL:KJZ97909.1}",1523156,Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;,N1_003_000G1_scaffold_0,Clostridia,Clostridiaceae,Clostridium,Clostridiales,Firmicutes,Clostridium sp. IBUN13A,Bacteria
N1_003_000G1_scaffold_0_11,A0A0F4WVU6_9CLOT,100.0,155,0,0,1,155,1,155,1.4999999999999996e-82,313.2,"Uncharacterized protein {ECO:0000313|EMBL:KJZ97869.1};  TaxID=1523156 species=""Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;"" source=""Clostridium sp. IBUN13A.;"" ",Uncharacterized protein {ECO:0000313|EMBL:KJZ97869.1},1523156,Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;,N1_003_000G1_scaffold_0,Clostridia,Clostridiaceae,Clostridium,Clostridiales,Firmicutes,Clostridium sp. IBUN13A,Bacteria
N1_003_000G1_scaffold_0_40,A0A0F4WWR7_9CLOT,100.0,208,0,0,1,208,1,208,1.2999999999999997e-113,416.8,"Phosphoserine phosphatase {ECO:0000313|EMBL:KJZ97910.1}; EC=3.1.3.3 {ECO:0000313|EMBL:KJZ97910.1};;  TaxID=1523156 species=""Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;"" source=""Clostridium sp. IBUN13A.;"" ",Phosphoserine phosphatase {ECO:0000313|EMBL:KJZ97910.1},1523156,Bacteria; Firmicutes; Clostridia; Clostridiales; Clostridiaceae; Clostridium.;,N1_003_000G1_scaffold_0,Clostridia,Clostridiaceae,Clostridium,Clostridiales,Firmicutes,Clostridium sp. IBUN13A,Bacteria
```

### Quick start:
```
./make_Tdb.py -h
usage: make_Tdb.py [-h] -b B6_LOC -o OUT_LOC [--version]

Make a Tdb.csv file

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

INPUT ARGUMENTS:
  -b B6_LOC, --b6_loc B6_LOC
                        location of b6+ file

OUTPUT ARGUMENTS:
  -o OUT_LOC, --out_loc OUT_LOC
                        location of output file
```

## tax_collector.py

This generates a taxonomy report from a b6+ file resulting from a single genome

### Quick start:
```
$ ./tax_collector.py -h
usage: tax_collector.py [-h] [-b B6_LOC] [-o OUT_LOC] [--SkipScaffolds]
                        [-stb SCAFFOLD2BIN] [--version]

Generate taxonomy report from a b6+ file resulting from a single genome

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

INPUT ARGUMENTS:
  -b B6_LOC, --b6_loc B6_LOC
                        location of b6+ file

OUTPUT ARGUMENTS:
  -o OUT_LOC, --out_loc OUT_LOC
                        output basename

OPTIONALL ARGUMENTS:
  --SkipScaffolds       skip generating per-scaffold taxonomy
  -stb SCAFFOLD2BIN, --scaffold2bin SCAFFOLD2BIN
                        scaffold to bin file for generating per-genome taxonomy
```

## quickTaxonomy_centrifuge.py

This program takes a genome bin, calls genes with prodigal, assigns those genes to an organism in a centrifuge database, and uses that information to determine the taxonomy of the bin.

The way that this functions is by mapping the nucleotide sequence of genes to a reference database. Thus, this doesn't work very well for novel organisms that aren't in the centrifuge database provided.

This program provides just a guess as to what the organism is, assuming that it is indeed in the reference database. That said, if an organism of the same species is in the centrifuge database, it works pretty well.

### Quick start:
```
$ ./quickTaxonomy_centrifuge.py -h
usage: quickTaxonomy_centrifuge.py [-h] [-ch CENT_HIT] [-cr CENT_REPORT]
                                   [-f FASTA] [--cent_index CENT_INDEX]
                                   [-m METHOD] [-per PERCENT]
                                   [--full_dump FULL_DUMP]
                                   [--min_score MIN_SCORE]
                                   [--min_diff MIN_DIFF] [--version]
                                   [-p PROCESSES]

Quickly determine the taxonomy of a genome

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -p PROCESSES, --processes PROCESSES
                        threads for centrifuge run

INPUT ARGUMENTS:
  -ch CENT_HIT, --cent_hit CENT_HIT
                        raw output of centrifuge run (hits.tsv)
  -cr CENT_REPORT, --cent_report CENT_REPORT
                        raw output of centrifuge run (report.tsv)
  -f FASTA, --fasta FASTA
                        location of fasta file(s) to detemine the taxonomy of
  --cent_index CENT_INDEX
                        path to centrifuge index (for example,
                        /home/mattolm/download/centrifuge/indices/b+h+v

OUTPUT ARGUMENTS:
  -m METHOD, --method METHOD
                        Method to determine taxonomy. percent = the lowest
                        taxonomic id with at least -p hits. max = the
                        taxonomic level with the most hits (lowest is genus
  -per PERCENT, --percent PERCENT
                        Minimum percent of genes for the percent method
  --full_dump FULL_DUMP
                        Output FULL taxonomy table to the file specified here
                        (for testing purposes)

FILTERING ARGUMENTS:
  --min_score MIN_SCORE
                        Minimum score for centrifuge hit
  --min_diff MIN_DIFF   Minimum score difference between first and second hit
```
