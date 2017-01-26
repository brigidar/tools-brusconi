#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	filter_reads.py								#
# Version     : 0.1									#
# Project     : SNPDV							#
# Description : Script to sort out no hits and get fasta		#
# Author      : Brigida Rusconi								#
# Date        : September 6th, 2016							#
#											#
#########################################################################################
#for replacement of a given value with NaN
#http://stackoverflow.com/questions/18172851/deleting-dataframe-row-in-pandas-based-on-column-value

# to remove any symbol and replace it with nan
#http://stackoverflow.com/questions/875968/how-to-remove-symbols-from-a-string-with-python

# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html

# for selecting rows that have an indel:
#http://stackoverflow.com/questions/14247586/python-pandas-how-to-select-rows-with-one-or-more-nulls-from-a-dataframe-without



#------------------------------------------------------------------------------------------


import argparse, os, sys, csv
import pdb
import numpy
from numpy import *
import Bio
from pandas import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="fasta snps", default="output_file.fasta")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-t', '--total', help="inverted table to look at", default= "snp_filtered_table.txt")
parser.add_argument('-r','--remove', help="remove non-canonical nucleotides", default= "False")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
output2_file = args.total
remove=args.remove

#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)
df=df.set_index(['molecule','refpos'])
#------------------------------------------------------------------------------------------
# count for each genome amount of no hits and report so that we know which genome contributed more to no hits
count_nohits=list()
for i,item in enumerate(df.columns):
    if 'qbase' in item:
        nohits=df.loc[:,item].value_counts()
        
        if 'No Hit' in nohits.index:
            count_nohits.append((item.split('qbase:')[1], str(nohits['No Hit'])))
        else:
            count_nohits.append((item.split('qbase:')[1], '0'))


with open('no_hits.txt','w') as output:
    for i in count_nohits:
        output.write('\t'.join(i) + '\n')




#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns
no_hit= df.mask(df=='No Hit')
removed=no_hit.dropna()
no_snp1=removed.mask(removed=='No SNP')
no_snp=no_snp1.dropna()
#pdb.set_trace()
print "No Hit removed: SNP left " + str(removed.index.size)
print "No SNP removed: SNP left " + str(no_snp.index.size)


bases=['A','C','G','T']
count_qbase=list(no_snp.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df2=no_snp.iloc[:,qindexes]
df4=no_snp.loc[:,'refbase']
df1=concat([df4,df2], axis=1, join_axes=[df2.index])
cols=df1.columns

# remove non-canononical snps (optional)
if remove=="True":
    df3=df1[~df1[cols].isin(bases).all(axis=1)].dropna(how='all')
    df1=df1[~df1.index.isin(df3.index)]
else:
    for i in bases:
        df1=df1[df1 !=i].dropna(how='all').fillna(i)
    print "identical lines removed: SNP left " + str( slash2.index.size)

removed1=df1.T
removed1.reset_index(inplace=True)
#-------------------------------------------------------------------------------------
#save file to fasta

tab2=[]
#pdb.set_trace()
for i in range(0,removed1.index.size):
    tab=[]

    if ':' in removed1.iloc[i,0]:
        tab.append(removed1.iloc[i,0].split(':')[1])
        tab.append(''.join(removed1.iloc[i,1:]))
    elif 'refbase' in removed1.iloc[i,0]:
        tab.append(removed1.iloc[i,0])
        tab.append(''.join(removed1.iloc[i,1:]))

    tab2.append(tab)

with open('table', 'w') as t:
    for i in range(0,len(tab2)):
        t.write("\t".join(tab2[i])+"\n")

with open('table','rU') as input:
    with open(output_file,'w') as output:
        sequences = SeqIO.parse(input, "tab")
        count = SeqIO.write(sequences, output, "fasta")
#pdb.set_trace()
#------------------------------------------------------------------------------------------
df=df[df.index.isin(df1.index)]
with open(output2_file,'w') as output2:
    df.to_csv(output2, sep='\t', index=False)

#------------------------------------------------------------------------------------------



