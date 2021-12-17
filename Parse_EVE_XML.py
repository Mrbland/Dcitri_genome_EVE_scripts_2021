"""
Script that removes doubles from a .xml BLASTx output and converts it to a csv file
Anneliek ter Horst & Curtis Carlson
17/12/21
"""


# imports
from __future__ import division
from Bio.Blast import NCBIXML
import csv
import sys
import os
from lxml import etree
import pandas as pd
from collections import defaultdict
from datetime import datetime



# Check if input file name was given
if len(sys.argv) > 2:
    result = NCBIXML.parse(open(sys.argv[1]))
    output = sys.argv[2]
else:
    print('Please specify the input file and make sure it is an xml file and also specify name for output files')
    raise SystemExit

# Write a header for the outputfile
header = ('sequence', 'length', 'perc_identity', 'gaps', 'frame', 'position_on_hit_start',
          'position_on_hit_stop', 'position_on_query_start', 'position_on_query_stop', 'evalue', 'score',  'query_id')

# open the outputfile
with open(output+'_raw.csv','w') as f:
    writer = csv.writer(f)
    writer.writerow(header)

    # Go into fasta records
    for record in result:

	    # Go into fasta alignments
	    if record.alignments:

	      # Check each alignment
	      for alignment in record.alignments:

	          # Make recognizable names for all xml input objects.
	          for hsp in alignment.hsps:
	            sequence = alignment.title
	            length = hsp.align_length
	            perc_identity = float((hsp.identities/hsp.align_length)*100)
	            gaps = hsp.gaps
	            query_frame = hsp.frame
	            query = record.query

	            # Hit is viral hit from viral database
	            position_on_hit_start = hsp.sbjct_start
	            position_on_hit_stop = hsp.sbjct_end

	            # Query is piRNA cluster of insect
	            position_on_query_start = hsp.query_start
	            position_on_query_stop = hsp.query_end
	            evalue = hsp.expect
	            score = hsp.score

	            # Write to csv
	            row = (sequence, length, perc_identity, gaps, query_frame[0],
	            position_on_hit_start, position_on_hit_stop ,position_on_query_start,
	            position_on_query_stop, evalue, score, query)
	            writer.writerow(row)

    # close the file
    f.close()
result.close()


# open the dataframe
df = pd.read_csv(open(sys.argv[2]+'_raw.csv'))
print('Number of entries from XML: ', len(df))


# max eval on position_on_query_start is equal
max_eval = df.groupby(['query_id', 'position_on_query_start']).score.transform(max)
df4 = df[df.score == max_eval]

# max eval on position_on_query_stop is equal
max_eval = df.groupby(['query_id', 'position_on_query_stop']).score.transform(max)
df5 = df[df.score == max_eval]

# merge both max tables
df = df4.append(df5)

# and remove where start sequence is equal
df = df.drop_duplicates(['query_id', 'position_on_query_start'])

# remove where stop sequence is equal
df = df.drop_duplicates(['query_id', 'position_on_query_stop'])

#remove where stop and start are equal
df = df.drop_duplicates([ 'query_id', 'position_on_query_start', 'position_on_query_stop'])

# output to csv
df.to_csv(sys.argv[2]+'_remove_dups.csv', index=False)

startTime = datetime.now()

## Functions
''' Function to iterate over the dataframe and remove the hits that are in each others frames, 
Keeping one with highest bitscore'''
def iterate_over_df(df):

    # list to store index that are either unique enough or have highest evalue
    results = []

    # list to save those that have already been added so they can be skiped
    to_be_skipped = []


    # loop through data frame
    for index, row in df.iterrows():

        # check if sequence or simmilar sequence already added
        if index in to_be_skipped:
            continue

        # initialize empty simmilar dict
        simmilar = defaultdict(int)

        results = unique(results)
        to_be_skipped = unique(to_be_skipped)

        for index2, row2 in df.iterrows():
        # check if position start or stop is equal and is not self.
            if index == index2:
                continue

            # check if possition start or stop is equal and is not self. (are they on same scaffold)
            #if row[11] == row2[11]:

            # if entry is comparing to itself
            if int(row[7]) == int(row2[7]) and int(row2[8]) == int(row[8]):
                continue

            elif (row[7] in range(int(row2[7]), int(row2[8])) or
            row2[7] in range(int(row[7]), int(row[8])) or
            row[8] in range(int(row2[7]), int(row2[8])) or
            row2[8] in range(int(row[7]), int(row[8]))):
                # add both indexes of simmilar sequences plus their score to the dict
                simmilar[index] = row[10]
                simmilar[index2] = row2[10]

        # check if simmilar sequences have been found
        if len(simmilar) > 0:

            # get the max evalue from the simmilar sequences
            max_index = max(simmilar.keys(), key=lambda k: simmilar[k])

            # add index with maximum evalue to results list
            results.append(max_index)

            # add checked indexs to be skipped list
            for selected_index, length in simmilar.items():
                    to_be_skipped.append(selected_index)

        # if seqeunce is unique add index to results
        if len(simmilar) == 0:
            results.append(index)
            to_be_skipped.append(index)
                
    df2 = df.loc[set(results)]
    return df2

#Function unique
'''Function to keep only unique units in a list. '''
def unique(seq):
    seen = set()
    return [seen.add(x) or x for x in seq if x not in seen]


# load dataframe
df = pd.read_csv(open(sys.argv[2]+'_remove_dups.csv'))
print('Number of entries after removing duplicates: ', len(df))

# Sort the df based on which scaffold it is
df = df.sort_values(by=['query_id'])

#create unique df of scaffolds
UniqueScaffolds = df.query_id.unique()

# create a list from that df
UniqueScaffolds_list = list(UniqueScaffolds)

#create a data frame dictionary to store your data frames
DataFrameDict = {elem : pd.DataFrame for elem in UniqueScaffolds}

for key in DataFrameDict.keys():
    DataFrameDict[key] = df[:][df.query_id == key]


# Make a list for all dfs that have more than 1 entry (otherwise they dont have to be checked)
unique_list_over_one =[]

# For each scaffold, check if more than one EVE is found
# If not, then we don't need to check for inframes. 
for scaffold in UniqueScaffolds_list:
    df = DataFrameDict[scaffold]
    if len(df) > 1:
        unique_list_over_one.append(scaffold)

# Check for scaffolds that have more than one EVE hit if there are any in frame
for scaffold in unique_list_over_one:
    DataFrameDict[scaffold] = iterate_over_df(DataFrameDict[scaffold])

# Concatenate all output dataframes into one
#df_all_concat = pd.concat([pd.concat(v) for k,v in DataFrameDict.items()])
df_all_concat = pd.concat(DataFrameDict ,axis=0, ignore_index=True)


print('Number of entries after removing in frame:', len(df_all_concat))

# Put the results in a new csv
df_all_concat.to_csv(sys.argv[2]+'_remove_inframe.csv', index=False)

print('Time taken: ', datetime.now() - startTime)


