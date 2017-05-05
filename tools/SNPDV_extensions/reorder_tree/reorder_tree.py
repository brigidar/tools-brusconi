#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	reorder_tree.py							#
# Version     : 1.1									#
# Project     : SNP Genotyper							#
# Description : Script to reorder SNP output by phylogeny.
# Author      : Brigida Rusconi								#
# Date        : March 10, 2016							#
#											#
#########################################################################################

# for isin information
#http://pandas.pydata.org/pandas-docs/stable/indexing.html




import argparse, os, sys, csv
#import pdb
from pandas import *
from Bio import Phylo


#output and input file name to give with the script and arguments for cutoff
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="SNPs reordered according to tree")
parser.add_argument('-s', '--snp_table', help="snp genotype output")
parser.add_argument('-t', '--tree', help="newick formatted tree")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
tree1 = args.tree

#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0)

#gets order of tree
def parse_tree(tree):
    names = []
    mytree = Phylo.read(tree, 'newick')
    tree_names = [ ]
    for clade in mytree.find_clades():
        if clade.name:
            names.append(clade.name)
    return names
#uses function parse_tree to get list of names
trees=parse_tree(tree1)

#get header with qbase
count_qbase=df.columns.values
qnames=[]
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qnames.append(v)
        qindexes.append(i)
# reorder snp file based on tree
df2= df.iloc[:,0:4]
for name in trees:
    for i in qnames:
        if i.split(':')[1]==name:
            df2=concat([df2,df.loc[:,i]], axis=1)
df3=concat([df2,df.iloc[:,(qindexes[-1]+1):]], axis=1)


#save reordered file
with open(output_file,'w') as output:
    df3.to_csv(output, sep='\t', index=False)






                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 


