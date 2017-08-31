#!/usr/bin/env python

#########################################################################################
#											#
# Name	      :	2_sort_2.py								#
# Version     : 0.6									#
# Project     : sort & merge SNP tables							#
# Description : Script to sort out no hits, indels, identical lines and double hits		#
# Author      : Brigida Rusconi								#
# Date        : August 4th, 2017							#
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
import argparse, os, sys, csv,pdb
#import pdb
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
parser.add_argument('-o', '--output', help="fasta snps")
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
df=df.set_index(['molecule','refpos']).fillna('--')


print "merged table number SNPS " + str(df.index.size)

#------------------------------------------------------------------------------------------
#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns
df=df.mask(df=='No Hit').dropna()

print "No Hit removed: SNP left " + str(df.index.size)

#------------------------------------------------------------------------------------------
# only columns with qbase and refbase in table
count_qbase=list(df.columns.values)

qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df1=df.iloc[:,qindexes]

#removes lines that had no hit because alignment was too short
df1=df1.mask(df1=='--').dropna()
ref=df['refbase']
df=df[df.index.isin(df1.index)]

#------------------------------------------------------------------------------------------
#multiple hits check for alignment length create dataframe with information

bindexes=[]
for i, v in enumerate(count_qbase):
    if 'blengths:' in v:
        bindexes.append(i)

df_l=df.iloc[:,bindexes]
maxim=[]
maxpos=[]

#split on / and find max value
for i in range(0, df_l.index.size):
    maxim2=[]
    for s in df_l.iloc[i,:]:
        if '/' in s:
            for n in s.split('/'):
                maxim2.append(int(n))
        else:
            maxim2.append(int(s))
    maxim.append(max(maxim2))

# get positions----------------------------------------------------

for i in range(0, df_l.index.size):
    maxpos2=[]
    for s in df_l.iloc[i,:]:
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

df2_r=df1.reset_index(drop=True, level=0)

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
            
df2_p=df2_r.set_index(df1.index)

#------------------------------------------------------------------------------------------
#creates dataframe with rows that have / in it by looking into each row

sl=[]
for i,v in enumerate(df2_p.index):
    for n in df2_p.iloc[i,:]:
        if len(n.split('/'))>1:
            sl.append(df2_p.index[i])

#------------------------------------------------------------------------------------------
# remove rows that are in empty

slash=df2_p.drop(sl)

print "double hits removed: SNP left " + str(slash.index.size)
slash2=concat([ref,slash], axis=1, join_axes=[slash.index])

#------------------------------------------------------------------------------------------
# remove identical line

bases=['A','C','G','T']
cols=slash2.columns

# remove non-canononical snps (optional)
if remove=="True":
    slash3=slash2[~slash2[cols].isin(bases).all(axis=1)].dropna(how='all')
    slash2=slash2[~slash2.index.isin(slash3.index)]
else:
    for i in bases:
        slash2=slash2[slash2 !=i].dropna(how='all').fillna(i)
    print "identical lines removed: SNP left " + str( slash2.index.size)

#------------------------------------------------------------------------------------------
#replaces lines with indel

indel=slash2.mask(slash2=='indel')
indel2=indel.dropna()
# selects the rows that have an indel to look at the frame shifts
frameshift=indel[indel.isnull().any(axis=1)]
frameshift2=df[df.index.isin(frameshift.index)]
#removes lines that have an indel and are intergenic
frameshift3=frameshift2.mask(frameshift2=='intergenic')
frameshift4=frameshift3.dropna()
print "Number of indel in genes " + str(frameshift4.index.size)
with open("frameshift_pos.txt",'w') as output:
    frameshift4.to_csv(output, sep='\t')

#------------------------------------------------------------------------------------------
# add transition and transversion information function developed by Mando Rodriguez

def get_trans_state(base, hit):
    
    if base == hit:
        return "--"
    elif base == 'A' and hit == 'G':
        return "transition"
    elif base == 'G' and hit == 'A':
        return "transition"
    elif base == 'C' and hit == 'T':
        return "transition"
    elif base == 'T' and hit == 'C':
        return "transition"
    
    elif base == 'A' and hit == 'C':
        return "transversion"
    elif base == 'C' and hit == 'A':
        return "transversion"
    
    elif base == 'A' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'A':
        return "transversion"
    
    elif base == 'C' and hit == 'G':
        return "transversion"
    elif base == 'G' and hit == 'C':
        return "transversion"
    
    elif base == 'G' and hit == 'T':
        return "transversion"
    elif base == 'T' and hit == 'G':
        return "transversion"
    
    else:
        return "Error [base:%s, hit:%s]" % (base, hit)
#------------------------------------------------------------------------------------------

indel3=indel2.mask(indel2=='NaN')
indel4=indel3.dropna()

#recalculate the transition transversions for new lists in case columns are removed.
#refbase list
df4=indel4['refbase']
ref_base=[]

for i in df4:
    ref_base.append(i)

#get alleles for each snp
qbas=indel4.iloc[:,1:]

snp_nb=[]
for i in range (0,qbas.index.size):
    snp_uniq=[n for n in unique(qbas.iloc[i,:])[:]]
    snp_nb.append(snp_uniq)

snp3=[]
for i,n in enumerate(snp_nb):
    snp4=[]
    for t in n:
        snp4+= t.split('/')
    snp3.append(snp4)

snp4=[]
for i,v in enumerate(snp3):
    snp2=[]
    for n in snp3[i]:
        if n!=ref_base[i]:
            snp2.append(n)
    snp4.append(snp2)

trs_trv=[]
for i,v in enumerate(ref_base):
    m=[n for n in snp4[i]]
    p=[]
    if len(m)>0:
        for t in m:
            p.append(get_trans_state(v,t))
        trs_trv.append('/'.join(p))
    else:
        trs_trv.append(get_trans_state(v,snp4[i][0]))

#------------------------------------------------------------------------------------------
print "removed lines with short alignment %s lines left" % (str(indel4.index.size))
final2 =indel4.reset_index(drop=True).T #drops indexes of molecule and refpos
final2=final2.reset_index()

#-------------------------------------------------------------------------------------
#save file to fasta

tab2=[]
for i in range(0,final2.index.size):
    tab=[]
    if ':' in final2.iloc[i,0]:
        tab.append(final2.iloc[i,0].split(':')[1])
        tab.append(''.join(final2.iloc[i,1:]))
    else:
        tab.append(final2.iloc[i,0])
        tab.append(''.join(final2.iloc[i,1:]))
    tab2.append(tab)

with open('table', 'w') as t:
    for i in range(0,len(tab2)):
        t.write("\t".join(tab2[i])+"\n")

with open('table','rU') as input:
    with open(output_file,'w') as output:
        sequences = SeqIO.parse(input, "tab")
        count = SeqIO.write(sequences, output, "fasta")

#------------------------------------------------------------------------------------------
#description rows back in overview table look at isin options with index it checks if the index is included in the other index and only keeps the one that are

fin=df[df.index.isin(indel4.index)]
fin2=fin.iloc[:,(qindexes[-1]+1):]
fin3=concat([indel4,fin2], axis=1, join_axes=[indel4.index])

#------------------------------------------------------------------------------------------
#condensate codon calling

codon=fin3['query_codon']

# create list with all the unique values for the codons
CC=[]
for i,v in enumerate(codon):
    if '/' in v:
        CC.append(v.split('/')[0])
    else:
        CC.append(v)

# SNP alleles are defined in list of lists snp4 and refbase in ref
# ref_codon column
ref_codon=[]
for i in fin3['ref_codon']:
    ref_codon.append(i)

#refbase list
ref_base=[]
for i in ref:
    ref_base.append(i)

#get position in codon
cod_pos=[]
for i,v in enumerate(ref_codon):
    if CC[i]=='--':
        cod_pos.append('--')
    else:
        if v[0]!=CC[i][0]:
            cod_pos.append(0)
        elif v[1]!=CC[i][1]:
            cod_pos.append(1)
        elif v[2]!= CC[i][2]:
            cod_pos.append(2)

#check if codon is still needed for the snp present
fin_codon=[]
for i,v in enumerate(ref_codon):
    if len(snp4[i])==1 and v!='--':
        t=list(ref_codon[i])
        t[cod_pos[i]]=snp4[i][0]
        fin_codon.append(''.join(t))
    elif v=='--':
        fin_codon.append('--')
    elif len(snp4[i])> 1:
        cd=[]
        for n in snp4[i]:
            t=list(v)
            t[cod_pos[i]]=n
            cd.append(''.join(t))
        fin_codon.append('/'.join(cd))

# with the set function it only keeps unique values because sometimes there is duplicates in the codons after merging

#------------------------------------------------------------------------------------------
#translate each codon to the corresponding aa using biopython

AA=[]
for i in fin_codon:
    codons = [n for n in i.split('/') if n != '--']
    # looks if there is only one empty one or an indel
    if len(codons) == 0:
        AA.append('--')
    else:
        aa = []
        #makes an amino acid out of each codon the if all() statement checks that all nucleotides are ATGC. If there is anything else the codon will not be translated and will be identified as unknown.
        for p in codons:
            if all([z in bases for z in p])==True:
                aa.append(str(Seq(p, generic_dna).translate(table=11)))
            else:
                aa.append('Unknown')
        AA.append('/'.join(aa))

#------------------------------------------------------------------------------------------
#replaces the query codons

for v,f in enumerate(fin3['query_codon']):
    fin3['query_codon'][v]=fin_codon[v]

#replaces the amino acids

for v,f in enumerate(fin3['query_aa']):
    fin3['query_aa'][v]=AA[v]

#replaces transitions/transversions:

for v,f in enumerate(fin3['transition/transversion']):
    fin3['transition/transversion'][v]=trs_trv[v]

#------------------------------------------------------------------------------------------
# compare query amino acids to reference amino acid to get syn to nsyn and intergenic

refaa=[]
fin4=fin3.reset_index()
for b in fin4['ref_aa']:
    refaa.append(b)

syn=[]
for i,v in enumerate(AA):
    if v=='--':
        syn.append('intergenic')
    elif '/' in v:
        amino = [n for n in v.split('/')]
        bb=[]
        for a in amino:
            if a==refaa[i]:
                bb.append('SYN')
            else:
                bb.append('NSYN')
        syn.append('/'.join(bb))
    else:
        if v==refaa[i]:
            syn.append('SYN')
        else:
            syn.append('NSYN')

for i,n in enumerate(syn):
    if n=='intergenic' and len(snp4[i])>1:
        syn[i]='/'.join(repeat(n,len(snp4[i])))

#insert column in dataframe
fin3.insert(0, 'syn?', syn)


fin3=fin3.reset_index(level=1)

#------------------------------------------------------------------------------------------
# put new values for snps_gene

snps_gene=fin3.groupby('gene_name').size().reset_index()

snps_gene2=[]
for i,v in enumerate(snps_gene.index):
    col=[n for n in snps_gene.iloc[i,:][:]]
    snps_gene2.append(col)

for i,v in enumerate(snps_gene2):
    for l,s in enumerate(fin3['gene_name']):
        if s=='intergenic':
            fin3['snps_per_gene'][l]='intergenic'
        elif s ==v[0]:
            fin3['snps_per_gene'][l]=snps_gene2[i][1]

fin3=fin3.reset_index()

for x,i in enumerate(fin3.snps_per_gene):
    if i != 'intergenic':
            fin3['snps/gene_length'][x]=float(i)/float(fin3.gene_length[x])

dn=fin3.groupby('gene_name')['syn?']

counter=list()
for name,group in dn:
    c1=0
    c2=0
    for x in group:
        if '/' in x:
            for n in x.split('/'):
                if n=='NSYN':
                    c1=c1+1
                else:
                    c2=c2+1
        else:
            if x=='NSYN':
                c1=c1+1
            else:
                c2=c2+1
    if c2 != 0:
        counter.append(repeat((c1/c2),len(group)).tolist())
    else:
        counter.append(repeat(0,len(group)).tolist())

flat_dn=[n for item in counter for n in item]
fin3['dn/ds'].update(Series(flat_dn))

for i,f in enumerate(fin3.columns):
    if 'num_hits' in f:
        fin3.drop(axis=1,inplace=True, labels=f)

#------------------------------------------------------------------------------------------
#save total file for plotting -t option
with open(output2_file,'w') as output2:
    fin3.to_csv(output2, sep='\t', index=False)








#for v,i in enumerate(snp4):
#    if i[0]=='A':
#        snp4[v][0]='T'
#    elif i[0]=='T':
#        snp4[v][0]='A'
#    elif i[0]=='C':
#        snp4[v][0]='G'
#    elif i[0]=='G':
#        snp4[v][0]='C'
#
#for v,i in enumerate(ref_base):
#    if i=='A':
#        ref_base[v]='T'
#    elif i=='T':
#        ref_base[v]='A'
#    elif i=='C':
#        ref_base[v]='G'
#    elif i=='G':
#        ref_base[v]='C'
#final=concat([df5,fin3],axis=1)
#final=concat([df5.T,fin3.T],axis=1)

