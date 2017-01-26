#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	drop_fasta.py								#
# Version     : 0.1									#
# Project     : drop a fasta						#
# Description : 		#
# Author      : Brigida Rusconi								#
# Date        : March 15th, 2016							#
#											#
#########################################################################################


import argparse, os, sys, csv
import pandas
#import pdb
import Bio
from pandas import *

from Bio import SeqIO

#------------------------------------------------------------------------------------------

#output and input file name to give with the script
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="filtered multifasta")
parser.add_argument('-f', '--multifasta', help="multifasta of snp")
parser.add_argument('-i', '--identifier',nargs='*', help="fasta identifier to remove")

args = parser.parse_args()
output_file = args.output
input_file = args.multifasta
x = args.identifier
y=[]
for n in x:
    y=n.split(',')

#------------------------------------------------------------------------------------------
excluded=[]
seq_iterator = []
for record in SeqIO.parse(open(input_file, 'rU'), "fasta"):
    if record.id in y:
        excluded.append(record)
    else:
        seq_iterator.append(record)

with open(output_file, 'w') as output_handle:
    SeqIO.write(seq_iterator, output_handle, "fasta")




