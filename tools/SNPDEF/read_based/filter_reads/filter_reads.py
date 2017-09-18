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

def get_snp(table):
    snp_n=[]
    ident=[]
    for i,item in enumerate(table.index):
        #only append snp
        snp_u=[n for n in unique(table.iloc[i,:])[:] if n!=table.refbase[i]]
        if len(snp_u)>0:
            snp_n.append(snp_u)
        else:
            snp_n.append([n for n in unique(table.iloc[i,:])[:]])
            ident.append(i)
    return snp_n, ident
#--------------------------------------------------------------------------------
def missing_char(str, pos,n):
    item=list(str)
    item[pos]=n
    str="".join(item)
    return str
#--------------------------------------------------------------------------------
# to invert nucleotide if gene is on opposite strand
def invert_nucl(nuc):
    nucleo=["A","G","T","C","N","R","Y", "S","W","K","M","B","V","D","H"]
    n_inv=["T","C","A","G","N","Y","R","S","W","M","K","V","B","H","D"]
    comb=zip(nucleo,n_inv)
    nucl_dict={}
    for nucleo,n_inv in comb:
        nucl_dict[nucleo]=n_inv
    nuc2=nucl_dict[nuc]
    return nuc2
#--------------------------------------------------------------------------------
bases=['A','C','G','T']
def get_trans_state(base, hit):
    if hit in bases:
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
        return "--"

#------------------------------------------------------------------------------------------



#output and input file name to give with the script
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="fasta snps", default="output_file.fasta")
parser.add_argument('-s', '--snp_table', help="snp table to sort")
parser.add_argument('-t', '--total', help="inverted table to look at", default= "snp_filtered_table.txt")
parser.add_argument('-r','--remove', help="remove non-canonical nucleotides", default= "False")
parser.add_argument('-c','--recalculate', help="Were any genomes removed?", default= "No")



args = parser.parse_args()
output_file = args.output
input_file = args.snp_table
output2_file = args.total
remove=args.remove
calc=args.recalculate
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

count_qbase=list(df.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)

print "Initial SNP count " + str(df.index.size)

#replaces lines with "No Hits" with NaN and removes lines with NaN in qbase columns
ex=df
ex=ex.replace({'No SNP':'Z'},regex=True)
ex=ex.replace({'No Hit':'Z'},regex=True)
exclude=ex[ex.apply(lambda row: row.astype(unicode).str.contains('Z', case=False).any(), axis=1)]
df=df[~df.index.isin(exclude.index)]

print "Missing SNP removed: SNP left " + str(df.index.size)


bases=['A','C','G','T']


df2=df.iloc[:,qindexes]
df4=df.loc[:,'refbase']
df1=concat([df4,df2], axis=1, join_axes=[df2.index])
cols=df1.columns

# remove non-canononical snps (optional)
if remove=="True":
    df3=df1[~df1[cols].isin(bases).all(axis=1)].dropna(how='all')
    df1=df1[~df1.index.isin(df3.index)]
    print "Non-canonical SNP removed: SNP left " + str( df1.index.size)
    df1=df1[df1 !=i].dropna(how='all').fillna(i)
    print "Identical positions removed: SNP left " + str( df1.index.size)
else:
    for i in bases:
        df1=df1[df1 !=i].dropna(how='all').fillna(i)
    print "Identical positions removed: SNP left " + str( df1.index.size)


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

#recalculate snps/gene gene length and dn/ds


#-------------------------------dn/ds-------------------------------
dn=df.groupby(['gene_name','syn?']).size().reset_index()
dn_2=dn[dn['syn?'].str.contains('SYN')]
dn_2.rename(columns={0:'count'},inplace=True)
dn_ds=dict()
    
for i,v in enumerate(dn_2['gene_name']):
    t=dn_2[dn_2['gene_name']==v]
    if any(t['syn?'].str.contains('/')):
        ns=dict()
        try:
            ns['NSYN']=t[t['syn?']=='NSYN']['count'].values[0]
        except IndexError:
            ns['NSYN']=0
        try:
            ns['SYN']=t[t['syn?']=='SYN']['count'].values[0]
        except IndexError:
            ns['SYN']=0
        bl=t[t['syn?'].str.contains('/')]
        if bl.index.size==1:
            c=bl.iloc[0,1].split('/')
            for d in c:
                if d=='NSYN':
                    ns['NSYN']=(ns['NSYN']+bl.iloc[0,2])
                else:
                    ns['SYN']=(ns['SYN']+bl.iloc[0,2])
        else:
            for g in range(0,bl.index.size):
                c=bl.iloc[g,1].split('/')
                for d in c:
                    if d=='NSYN':
                        ns['NSYN']=(ns['NSYN']+bl.iloc[g,2])
                    else:
                        ns['SYN']=(ns['SYN']+bl.iloc[g,2])
        dn_ds[v]=ns['NSYN']/ns['SYN']
    else:
        try:
            dn_ds[v]=t[t['syn?']=='NSYN']['count'].values[0] / t[t['syn?']=='SYN']['count'].values[0]
        except IndexError:
            dn_ds[v]=0

df.drop(['snps_per_gene','dn_ds','snps/gene_length','strand'],axis=1,inplace=True)
df['dn_ds']=df['gene_name'].map(dn_ds)


snps_gene=df.groupby('gene_name').size().reset_index()
snps_gene2=dict(zip(snps_gene['gene_name'].tolist(),snps_gene[0].tolist()))
df['snps_per_gene']=df['gene_name'].map(snps_gene2)
df['snps/gene_length']=df['snps_per_gene'].astype(float, errors='ignore')/df['gene_length'].astype(float,errors='ignore')

if calc=='No':
    with open(output2_file,'w') as output2:
        df.to_csv(output2, sep='\t')


#-------------------------------------------------------------------

else:

    cod=df.dropna(subset=['gene_name'])
    cod.reset_index(inplace=True)

    count_qbase2=list(cod.columns.values)
    qindexes2=[]
    for i, v in enumerate(count_qbase2):
        if 'qbase:' in v:
            qindexes2.append(i)

    # get query base information
    df3=cod.iloc[:,qindexes2].join(cod.refbase)
    #position in codon
    pos1=(cod.pos_in_gene.astype(int) % 3).tolist()
    pos1=[ x if x!=0 else 3 for x in pos1 ]
    pos1=[(x-1) for x in pos1]
    ref_codon=cod.ref_codon.astype(str).tolist()
    #get allele for each position


    snp_nb, idn =get_snp(df3)


    query_codon=[]
    for i,v in enumerate(snp_nb):
        if ref_codon[i]=='nan':
            query_codon.append('nan')
        else:
            #positions on -1 strand need to be inverted use invert_nucl function
            if cod['strand'][i]==1:
                ts=list()
                if len(v)==1:
                    query_codon.append(missing_char(str(ref_codon[i]),pos1[i],v[0]))
                else:
                    #multiallelic position gets codon for each
                    for n in v:
                        ts.append(missing_char(str(ref_codon[i]),pos1[i],n[0]))
                    query_codon.append('/'.join(ts))
            else:
                ts=list()
                if len(v)==1:
                    query_codon.append(missing_char(str(ref_codon[i]),pos1[i],invert_nucl(v[0])))
                else:
                            #multiallelic position gets codon for each
                    for n in v:
                        ts.append(missing_char(str(ref_codon[i]),pos1[i],invert_nucl(n[0])))
                        query_codon.append('/'.join(ts))
    query_aa=[]
    for i,v in enumerate(query_codon):
        qq=[]
        if v=='nan':
            query_aa.append('nan')
        else:
            if '/' in v:
                gl=v.split('/')
                for n in gl:
                    qq.append(str(Seq(n, generic_dna).translate(table=t_t)))
                query_aa.append('/'.join(qq))
            else:
                query_aa.append(str(Seq(v, generic_dna).translate(table=t_t)))

    cod.insert(cod.columns.size,'query_codon',query_codon)
    cod.insert(cod.columns.size,'query_aa',query_aa)

    print "Read query codons and aa"
    #-------------------------------transition/transversion -------------------------------
    #
    df4=df.iloc[:,qindexes].join(df.refbase)


    snp_nb2, ident=get_snp(df4)
    ts_tv=[]
    for i,v in enumerate(snp_nb2):
        ts=[]
        if len(v)==1 and v[0]!= 'No Hit':
            ts_tv.append(get_trans_state(df.refbase[i],v[0]))
        else:
            for n in v:
                if n!='No Hit':
                    ts.append(get_trans_state(df.refbase[i],n))
            ts_tv.append('/'.join(ts))
    df.insert(df.columns.size,'transition/transversion',ts_tv)
    print "Read transition/transversion"
    df.drop('refpos_norm',axis=1,inplace=True)


    cod.set_index(['molecule','refpos'],inplace=True)
    fin=df.join(cod)
    # -------------------------------synonymous nonsynonymous-------------------------------
    query_aa=fin.query_aa.astype(str).tolist()
    ref_aa=fin.ref_aa.astype(str).tolist()
    syn=[]
    for i,item in enumerate(query_aa):
        if item=='nan':
            if len(snp_nb2[i])==1:
                syn.append('intergenic')
            else:
                g=[n for n in snp_nb2[i] if n!='No Hit']
                syn.append('/'.join(repeat('intergenic',len(g))))
        elif '/' in item:
            mult=list()
            for n in item.split('/'):
                if n==ref_aa[i]:
                    mult.append('SYN')
                else:
                    mult.append('NSYN')
            syn.append('/'.join(mult))
        else:
            if item==ref_aa[i]:
                syn.append('SYN')
            else:
                syn.append('NSYN')

    fin.insert(0,'syn?',syn)

    #genes that are not CDS are tagged as genic
    genic=dict()
    genic['No CDS']='genic'
    fin.reset_index(inplace=True)
    fin.set_index('product',inplace=True)
    fin.replace({'syn?':genic},inplace=True)
    fin.reset_index(inplace=True)

    #SNPs that are actually identical are replaced with No SNP
    fin['syn?'][ident]='No SNP'

    #-------------------------------dn/ds-------------------------------
    dn=fin.groupby(['gene_name','syn?']).size().reset_index()
    dn_2=dn[dn['syn?'].str.contains('SYN')]
    dn_2.rename(columns={0:'count'},inplace=True)
    dn_ds=dict()

    for i,v in enumerate(dn_2['gene_name']):
        t=dn_2[dn_2['gene_name']==v]
        if any(t['syn?'].str.contains('/')):
            ns=dict()
            try:
                ns['NSYN']=t[t['syn?']=='NSYN']['count'].values[0]
            except IndexError:
                ns['NSYN']=0
            try:
                ns['SYN']=t[t['syn?']=='SYN']['count'].values[0]
            except IndexError:
                ns['SYN']=0
            bl=t[t['syn?'].str.contains('/')]
            if bl.index.size==1:
                c=bl.iloc[0,1].split('/')
                for d in c:
                    if d=='NSYN':
                        ns['NSYN']=(ns['NSYN']+bl.iloc[0,2])
                    else:
                        ns['SYN']=(ns['SYN']+bl.iloc[0,2])
            else:
                for g in range(0,bl.index.size):
                    c=bl.iloc[g,1].split('/')
                    for d in c:
                        if d=='NSYN':
                            ns['NSYN']=(ns['NSYN']+bl.iloc[g,2])
                        else:
                            ns['SYN']=(ns['SYN']+bl.iloc[g,2])
            dn_ds[v]=ns['NSYN']/ns['SYN']
        else:
            try:
                dn_ds[v]=t[t['syn?']=='NSYN']['count'].values[0] / t[t['syn?']=='SYN']['count'].values[0]
            except IndexError:
                dn_ds[v]=0

    fin['dn_ds']=fin['gene_name'].map(dn_ds)
    fin1=fin.iloc[:,1:(max(qindexes)+3)]
    fin1.set_index(['molecule','refpos'],inplace=True)
    fin2=fin.reindex_axis(['molecule','refpos','gene_name','gene_start','gene_end','gene_length','pos_in_gene','ref_codon','ref_aa','query_codon','query_aa','product','transition/transversion','snps_per_gene','snps/gene_length','dn_ds','strand'],axis=1)#'dn/ds'
    fin2.set_index(['molecule','refpos'],inplace=True)
    final=fin1.join(fin2)
    final.reset_index(inplace=True)
    final.sort_values(by=['molecule','refpos'],inplace=True)
    with open(output2_file,'w') as output2:
        final.to_csv(output2, sep='\t',index=False)

