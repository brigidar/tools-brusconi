#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	plots_2.py								#
# Version     : 0.3									#
# Project     : plot snp according to sliding window							#
# Description : provide size of genome and sliding window to plot the snps along the chromosome axis		#
# Author      : Brigida Rusconi								#
# Date        : March 17th, 2016							#
#											#
#########################################################################################
# sliding window function definition from old python docs
#http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python

import argparse, os, sys, csv
#import pdb
from pandas import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import islice
#length of chromosome to have a x scale and input file
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--length_chr', help="length of chromosome")
parser.add_argument('-s', '--snp_table', help="snp table for snp distribution")
parser.add_argument('-b', '--bin', help="bin number for barplot")

args = parser.parse_args()
input_file = args.snp_table
chr=args.length_chr
num_bins=args.bin

# read in a genbank file and make a dictionary of length and molecule name if you have multiple contigs

chr=int(chr)
bins=int(num_bins)


#sliding window function
def window2(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    sl=[]
    result = tuple(islice(it, n))
    if len(result) == n:
        yield sum(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield sum(result)


#function that chops list in equal pieces
#def chunks(positions, bin):
#for i in xrange(0, len(positions), bin):
#  yield positions[i:i+bin]

#read in file as dataframe
df =read_csv(input_file, sep='\t', dtype=unicode)

#check in each position if there is a SNP if yes convert from 0 to 1 in the pos list
for name,group in df.groupby('molecule'):
    pos=[0 for i in range(0,chr)]
    df2=DataFrame(range(0,chr), columns=['refpos'])
    chunk=[]
    interval=[]
    for i,t in enumerate(group['refpos']):
        if int(t) in range(0,chr):
            pos[int(t)]=1
    df6=DataFrame(pos, columns=['snp_pos'])
    df7=concat([df2,df6], axis=1, join_axes=[df2.index])
    sl=[]
#the number of snps is counted for all positions within the sliding window size
    for each in window2(pos,n=bins):
        sl.append(each)


    interval2=range(bins,(len(sl)+bins))

    df3=DataFrame(interval2, columns=['refpos'])
    df3.insert(df3.columns.size,'snp_window',sl)
    with open('plots.txt', 'w') as table:
        df3.to_csv(table, sep='\t')
    fig=df3.plot(kind='area', x='refpos',y='snp_window', title=name, legend=False, figsize=(15,5))
    fig2=fig.get_figure()
    fig2.savefig('plot.png')





#plt.show()




