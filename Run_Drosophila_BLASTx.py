"""
Script that removes hits to the Drosophila melanogaster proteome and outputs csv 
Anneliek ter Horst & Curtis Carlson
17/12/21
"""


# load modules
import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
from datetime import datetime
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastxCommandline

# Check if input file name was given
if len(sys.argv) > 1:
    insect = sys.argv[1]
    df = pd.read_csv(open('../'+sys.argv[1]+'/'+sys.argv[1]+'.csv'))
    print(len(df))
else:
    print('Please specify the input file and make sure it is a csv file')
    raise SystemExit


# set current directory, where scripts are
scripts_dir = Path.cwd()

# Set parent directory
parent_dir = str(Path().resolve().parent)

# Make new dir with bed files
working_dir = parent_dir+'/'+insect

# Now navigate to working dir
os.chdir(working_dir)
print(working_dir)

# Make all the folders needed later

# Make a directory for the bed file
print('Making directories')
if not os.path.exists('bed_files'):
    os.mkdir('bed_files')

# Make a directory for nucleotide sequences
if not os.path.exists('nt_sequences'):
    os.mkdir('nt_sequences')

# Make directory for blast results against drosophila proteome
if not os.path.exists('BLAST_results'):
    os.mkdir('BLAST_results')


# Now for the bed files, we need:
# contig name (aka direction), pos on query start, and stop
df_bed_file = df[['position_on_query_start', 'position_on_query_stop', 'query_id']]

# Make the floats into ints, so that bedtools wont complain
cols = ['position_on_query_start', 'position_on_query_stop']
df_bed_file[cols] = df_bed_file[cols].applymap(np.int64)

# Biologists use spaces in their names so replace the spaces with underscores, otherwise youll have problems,
df_bed_file = df_bed_file.replace(' ', '_', regex=True)

# save as tsv (tab separated), no header and no index column cause that screws bed files
# Just saving this whole file so that we dont have to redo all this
df_bed_file.to_csv('./bed_files/'+insect+'.bed', sep="\t", index=False, header=False)
print('Made bed file')
# read names and postions from bed file
positions = defaultdict(list)
with open('./bed_files/'+insect+'.bed') as f:
    for line in f:
        start, stop, name = line.split("\t")
        name = name.strip('\n')
        positions[name].append((int(start), int(stop)))

# parse faste file and turn into dictionary
records = SeqIO.to_dict(SeqIO.parse(open('./genome/'+insect+".fa"), 'fasta'))
records
# search for short sequences
short_seq_records = []
for name in positions:
    for (start, stop) in positions[name]:
        long_seq_record = records[name]
        long_seq = long_seq_record.seq
        #alphabet = long_seq.alphabet
        short_seq = str(long_seq)[start-1:stop]
        short_seq_record = SeqRecord(Seq(short_seq), id=name+'.'+str(start)+'.'+str(stop), description='')
        short_seq_records.append(short_seq_record)

# # write to file
with open('./nt_sequences/'+insect+'.fasta', 'w') as f:
    SeqIO.write(short_seq_records, f, 'fasta')
print('made fasta file')
# Do blast with biopython
# Define the blast command
print('now running drosophila proteome blast')

blastx_cline = NcbiblastxCommandline(query='./nt_sequences/'+insect+'.fasta',
                                     db='./BLAST_db/drosophila_no_spaces.fasta',
                                     evalue=0.001, outfmt=10, max_hsps=1, num_threads=6,
                                     max_target_seqs=1, out="./BLAST_results/"+insect+'.csv')
# Run blast
stdout, stderr = blastx_cline()


# Do something with pandas to concatenate and keep
# the ones that do not have a match with drosophila proteome
drosophila_blast_df = pd.read_csv('./BLAST_results/'+insect+'.csv',
                  sep=',',
                  names=['query_id', 'subject_id', 'perc_identity','length', 'mismatches', 'gap_open', 'position_on_hit_start',
                         'position_on_hit_stop', 'position_on_query_start',
                         'position_on_query_stop','evalue',
                         'score',  'direction','x','y'])

print('done with drosophila proteome blast')


drosophila_blast_df = drosophila_blast_df[['query_id']]

drosophila_blast_df[['query_id','position_on_query_start','position_on_query_stop']] = drosophila_blast_df['query_id'].str.split('.',expand=True)

drosophila_blast_df.to_csv('./BLAST_results/'+insect+'_drosophila_hits.csv', index=False)


# merge the drosophila output and original df output
df['position_on_query_start']=df['position_on_query_start'].astype(int)
df['position_on_query_stop']=df['position_on_query_stop'].astype(int)

drosophila_blast_df['position_on_query_start']=drosophila_blast_df['position_on_query_start'].astype(int)
drosophila_blast_df['position_on_query_stop']=drosophila_blast_df['position_on_query_stop'].astype(int)

# Merge the drosophila blast file and the viral blast hits on the query, start and stop
common = df.merge(drosophila_blast_df,on=['query_id', 'position_on_query_start','position_on_query_stop'])
print(len(common))

# Keep only hits in viral blast that do not have a hit in drosophila proteome
df_merged = pd.concat([df, common, common]).drop_duplicates(keep=False)

print(len(df))
print(len(df_merged))
df_merged.to_csv(insect+'_remove_drosophila_hits.csv', index=False)
