#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	filter_nucmer.py								#
# Version     : 0.1									#
# Project     : sort & merge SNP tables							#
# Description : Script to extract nucmer repeats		#
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
#import pdb
import numpy as np
from pandas import *

#------------------------------------------------------------------------------------------


#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="sorted nucmer repeats")
parser.add_argument('-n', '--nucmer', help="nucmer")


args = parser.parse_args()
output_file = args.output
input_file = args.nucmer
#------------------------------------------------------------------------------------------


#read in file as dataframe
df =read_csv(input_file,sep='\t', dtype=object)

#------------------------------------------------------------------------------------------
for i in range(0,df.index.size):
    if df.iloc[i,0]==df.iloc[i,2] and df.iloc[i,1]==df.iloc[i,3]:
        df.iloc[i,0]=0

df2=df[df.iloc[:,0]>0]
#pdb.set_trace()
#repeat regions
df3=concat([df2.iloc[:,7],df2.iloc[:,0:2]], axis=1)
df3=df3.T.reset_index(drop=True).T
#repeat regions check orientation
df4=concat([df2.iloc[:,8],df2.iloc[:,2:4]], axis=1)

#check if needs to be inverted
new_index=[]
new_index2=[]
for i in range(0,df4.index.size):
    if df4.iloc[i,1]>df4.iloc[i,2]:
        new_index.append(df4.index[i])
    else:
        new_index2.append(df4.index[i])

#check if there are any inverted repeats
if len(new_index)>0:
    df5=df4[df4.index.isin(new_index)]
    #pdb.set_trace()
    cols = df5.columns.tolist()
    #pdb.set_trace()
    cols=[cols[0]]+cols[-1:]+[cols[1]]
    df5=df5[cols].T.reset_index(drop=True).T
    #normal order
    if len(new_index2)>0:
        df6=df4[df4.index.isin(new_index2)].T.reset_index(drop=True).T
        df7=df5.append(df6)
        df8=df7.append(df3)
    #inverted order fixed
    else:
        df8=concat([df5,df3])
#pdb.set_trace()
else:
    df4=df4.T.reset_index(drop=True).T
    df8=concat([df3,df4])
#------------------------------------------------------------------------------------------



#save nucmer repeats for exclusion
with open(output_file,'w') as output2:
    df8.to_csv(output2, sep='\t', index=False, header=None)


