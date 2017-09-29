# tRep
Quickly get the taxonomy of a genome.

### OVERVIEW
The program takes a genome bin, calls genes with prodigal, assigns those genes to an organism in a centrifuge database, and uses this ifnormation to determine the taxonomy of the bin.

### WARNING
The way that this functions is by mapping the nucleotide sequence of genes to a reference database. Thus, this doesn't work very well for novel organisms that aren't in the centrifuge database provided. 

This program provides just a guess as to what the organism is, assuming that it is indeed in the reference database. That said, if an organism of the same species is in the centrifuge database, it works pretty well.

### Dependencies

dRep (https://github.com/MrOlm/drep)

ETE3 (https://github.com/etetoolkit/ete)

### Quick start:
```
$ ./tRep.py -h
usage: tRep.py [-h] [-ch CENT_HIT] [-cr CENT_REPORT] [-f FASTA]
               [--cent_index CENT_INDEX] [-m METHOD] [-per PERCENT]
               [--full_dump FULL_DUMP] [--min_score MIN_SCORE]
               [--min_diff MIN_DIFF] [--version] [-p PROCESSES]

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
