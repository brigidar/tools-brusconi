#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	filter_blastn.py								#
# Version     : 0.1									#
# Project     : sort & merge SNP tables							#
# Description : find blsatn hits of IS elements in reference genome	#
# Author      : Brigida Rusconi								#
# Date        : March 9th, 2016							#
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
import numpy as np
from pandas import *
#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="sorted filtered table for fasta", default="output.txt")
parser.add_argument('-b', '--blastn', help="blastn input")


args = parser.parse_args()
output_file = args.output
input_file = args.blastn
#------------------------------------------------------------------------------------------
#read in file as dataframe

#checks if there are any blast results before continuing
#https://stackoverflow.com/questions/43197556/parsing-empty-file-with-no-columns

from pandas.io.common import EmptyDataError

try:
    df =read_csv(input_file,sep='\t')
except EmptyDataError:
    df = DataFrame()

if df.empty==True:
    with open(output_file,'w') as output2:
        df.to_csv(output2)

else:
    for i in range(0,df.index.size):
        if float(df.iloc[i,3])/float(df.iloc[i,23])<0.5:
            df.iloc[i,3]=0

    df2=df[df.iloc[:,3]>0]

    df3=concat([df2.iloc[:,0],df2.iloc[:,6:8]], axis=1)
    new_index=[]
    new_index2=[]
    for i in range(0,df3.index.size):
        if df3.iloc[i,1]>df3.iloc[i,2]:
            new_index.append(df3.index[i])
        else:
            new_index2.append(df3.index[i])


    df5=df3[df3.index.isin(new_index2)].T.reset_index(drop=True).T
    #inverted matches (start bigger than stop)
    if len(new_index)>0:
        df4=df3[df3.index.isin(new_index)]

        cols = df4.columns.tolist()

        cols=[cols[0]]+cols[-1:]+[cols[1]]
        df4=df4[cols].T.reset_index(drop=True).T
    #normal order


        df6=df5.append(df4)

    else:
        df6=df5

    #------------------------------------------------------------------------------------------




    #save IS regions bed format
    with open(output_file,'w') as output2:
        df6.to_csv(output2, sep='\t', index=False, header=None)


