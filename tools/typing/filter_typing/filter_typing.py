#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	filter_typing.py								#
# Version     : 0.1									#
# Project     : sort & merge SNP tables							#
# Description : Script to extract multifasta from typing output		#
# Author      : Brigida Rusconi								#
# Date        : March 15th, 2016							#
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
#import pdb
import numpy as np
from pandas import *
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="fasta snps")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-p', '--positions', help="positions")

args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
input_file2=args.positions
#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)
#need to fill na otherwise it cannot do boolean operations
#pdb.set_trace()
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------

# only columns with qbase and refbase in table
count_qbase=list(df.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df2=df.iloc[:,qindexes]

#replaces no hit and possible indels as N
for i in range(0,df2.index.size):
    for n in df2.iloc[i,:]:
        if n=="No Hit" or n=="indel":
            n="N"

#------------------------------------------------------------------------------------------

#mutiple hits check for alignment length to select if keeping or not
bindexes=[]
for i, v in enumerate(count_qbase):
    if 'blengths:' in v:
        bindexes.append(i)
df_l=df[df.index.isin(df2.index)]
df_l2=df_l.iloc[:,bindexes]

maxim=[]
maxpos=[]
#split on / and find max value
for i in range(0, df_l2.index.size):
    maxim2=[]
    for s in df_l2.iloc[i,:]:
        
        if '/' in s:
            for n in s.split('/'):
                maxim2.append(int(n))
        else:
            maxim2.append(int(s))
    maxim.append(max(maxim2))


# get positions----------------------------------------------------

for i in range(0, df_l2.index.size):
    maxpos2=[]
    for s in df_l2.iloc[i,:]:
        if '/' in s:
            maxpos3=[]
            for n,m in enumerate(s.split('/')):
                if int(m)==maxim[i]:
                    maxpos3.append(n)
            maxpos2.append(maxpos3)
        

        else:
            maxpos3=['0']
            maxpos2.append(maxpos3)
    maxpos.append(maxpos2)

# check back positions with nucleotides----------------------------------------------
df2_r=df2.reset_index(drop=True, level=0)

for i,v in enumerate(maxpos):
    
    for a,b in enumerate(v):
        if '--' in b:
            df2_r.iloc[i,a]='--'
    
        if '/' in df2_r.iloc[i,a]:
            g=[]
            m=[]
            m=df2_r.iloc[i,a].split('/')
            if len(b)==1:
                df2_r.iloc[i,a]=m[int(b[0])]
            
            elif len(b)==0:
                df2_r.iloc[i,a]='NaN'

            else:
                
                for c in b:
                    g.append(m[c])
                df2_r.iloc[i,a]=str('/'.join(g))
            
df2_p=df2_r.set_index(df2.index)




#------------------------------------------------------------------------------------------
sl=[]
#creates dataframe with rows that have / in it by looking into each row and if the item has a / it will
for i,v in enumerate(df2_p.index):
    for n in df2_p.iloc[i,:]:
        if len(n.split('/'))>1:
            
            sl.append(df2_p.index[i])


#------------------------------------------------------------------------------------------

# remove rows that are in empty
slash=df2_p.drop(sl)

df3=read_csv(input_file2,sep='\t', dtype=object, names=['molecule','refpos'])
df4=df.reindex(slash.index)

df5=concat([df4.iloc[:,0:4], slash], axis=1)
df6=concat([df5, df4.iloc[:, qindexes[-1]+1:]], axis=1)
df6=df6.set_index(df6['refpos'])
df7=df6.reindex(df3['refpos'])
df7=df7.fillna('N')
with open("typing.txt",'w') as output:
    df7.to_csv(output, sep='\t', index=False)

#-------------------------------------------------------------------------------------
#save file to fasta
df7=df7.reset_index(drop=True)
df8=df7.iloc[:,qindexes].T.reset_index()

tab2=[]
#pdb.set_trace()
for i in range(0,df8.index.size):
    tab=[]

    if ':' in df8.iloc[i,0]:

        tab.append(df8.iloc[i,0].split(':')[1])
        tab.append(''.join(df8.iloc[i,1:]))
    else:
        tab.append(df8.iloc[i,0])
        tab.append(''.join(df8.iloc[i,1:]))


    tab2.append(tab)
with open('table', 'w') as t:
    for i in range(0,len(tab2)):
        t.write("\t".join(tab2[i])+"\n")

with open('table','rU') as input:
    with open(output_file,'w') as output:
        sequences = SeqIO.parse(input, "tab")
        count = SeqIO.write(sequences, output, "fasta")

#------------------------------------------------------------------------------------------





