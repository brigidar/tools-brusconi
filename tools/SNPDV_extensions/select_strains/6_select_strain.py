#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	6_select_strain.py								#
# Version     : 1.12									#
# Project     : SNP tables downstream analysis							#
# Description : this script allows you to give a list of strain names and get the genotype profiles that match that list		#
# Author      : Brigida Rusconi								#
# Date        : March 10th, 2016							#
#											#
##########################################################################################

import argparse, os, sys, csv
#import pdb
from pandas import *
#------------------------------------------------------------------------------------------

#output and input file name
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="groups that contain that genome")
parser.add_argument('-s', '--genotyper_table', help="genotyper summary file", default="summary_groups.txt")
parser.add_argument('-g', '--genome',nargs='*', help="strain name list separated by commas")


args = parser.parse_args()
output_file = args.output
input_file = args.genotyper_table
x = args.genome
y=[]
for n in x:
    y=n.split(',')

#------------------------------------------------------------------------------------------

#read in files as dataframe
df =read_csv(input_file,header=0, sep='\t', dtype=unicode)

#------------------------------------------------------------------------------------------

# select rows that have the genome of interest in there
df1=df.copy()
for i in y:
    df1=df1[df1['Genomes_w_1'].str.contains(str(i))]

# limit to rows that have only the list items and nothing more

df1b=DataFrame()
for b in df1.index:
    if len(df1.loc[b,'Genomes_w_1'].split(','))==len(y):
        df1b=concat([df1b,df1.loc[b,:]], axis=1)

df1b=df1b.T

#------------------------------------------------------------------------------------------


df2=df.copy()
for i in y:
    df2=df2[df2['Genomes_w_2'].str.contains(str(i))]

df2b=DataFrame()
for b in df2.index:
    if len(df2.loc[b,'Genomes_w_2'].split(','))==len(y):
        df2b=concat([df2b,df2.loc[b,:]], axis=1)
df2b=df2b.T

#------------------------------------------------------------------------------------------

df3=df.copy()
for i in y:
    df3=df3[df3['Genomes_w_3'].str.contains(str(i))]

df3b=DataFrame()
for b in df3.index:
    if len(df3.loc[b,'Genomes_w_3'].split(','))==len(y):
        df3b=concat([df3b,df3.loc[b,:]], axis=1)
df3b=df3b.T

#------------------------------------------------------------------------------------------

df4=df.copy()

for i in y:
    df4=df4[df4['Genomes_w_zero'].str.contains(str(i))]
df4b=DataFrame()
for b in df4.index:
    if len(df4.loc[b,'Genomes_w_zero'].split(','))==len(y)+1:
        df4b=concat([df4b,df4.loc[b,:]], axis=1)
df4b=df4b.T
#pdb.set_trace()
#------------------------------------------------------------------------------------------

df5=concat([df1b,df2b], axis=0)
df6=concat([df4b, df3b], axis=0)
df7=concat([df5, df6], axis=0)
#------------------------------------------------------------------------------------------
#pdb.set_trace()
#------------------------------------------------------------------------------------------


#save selected table
with open(output_file ,'w') as output:
    if df7.index.size==0:
        output.write( "No genotype group for this list \n")
    else:
        df7.to_csv(output, sep='\t')



