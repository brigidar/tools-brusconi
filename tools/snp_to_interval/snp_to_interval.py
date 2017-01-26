#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	snp_to_interval.py								#
# Version     : 0.1								#
# Project     : make interval from SNPs							#
# Description : make interval from SNPs to exclude		#
# Author      : Brigida Rusconi								#
# Date        : March 11th, 2016							#
#											#
#########################################################################################


import argparse, os, sys, csv
#import pdb
from pandas import *



#output file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-s', '--input', help="SNP predicted")
parser.add_argument('-o', '--output', help="output file", default="output.txt")

args = parser.parse_args()
input_file = args.input
output=args.output

df =read_csv(input_file,sep='\t', header=None, dtype=object)

# skips first column with chr or locus_tag

table2=[]
for i in range(0,df.index.size):
    table2.append(int(df.iloc[i,1])+1)

df.insert(2,'stop',table2)



with open(output,'w') as of:
     df.to_csv(output, sep='\t', index=False, header=None)