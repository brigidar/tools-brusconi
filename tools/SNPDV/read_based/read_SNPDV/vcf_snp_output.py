#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	vcf_snp_output.py								#
# Version     : 0.3								#
# Project     : snp verify reads						#
# Description : Script to populate SNPs identified from reads with location information		#
# Author      : Brigida Rusconi								#
# Date        : January 9th, 2017						#
#											#
#########################################################################################


import argparse, os, sys, csv
import pdb
#import glob,logging
from pandas import *
from numpy import *
#from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SearchIO
from Bio.SeqFeature import FeatureLocation
#------------------------------------------------------------------------------------------
# got the initial functions from Mando's script had to adapt the rest as I already have a table with the SNPs.
# There is only one Genbank file read in as there is only one reference.
#-----------------------------------------------------------------------

# Description   : Parse reference GenBank files to obtain the gene information for reference genome
# Parameters    : gbk_files = Path to the list of files of reference GenBank genomes
# Returns       : retval = Reference to a hash of gene information. Key = Molecule name, gene name, tags like start, end, length, strand, gene, pseudo, product
#	Value = Annotation of each gene
#	seq_cache = Reference to hash of reference genome sequences. Key = Molecule name (display id from GenBank file)
#	Value = FASTA sequence

#------------------------------Functions-------------------------------------------------

def parse_genbank_file_list(genbank_file):
    
    genbank_recs = []
        
    genbank_input_handle = open(genbank_file, "rU")
    for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):
        
        genbank_recs.append(gbrec)

    genbank_input_handle.close()
    
    print "Read in %i genbank records from file %s" % (len(genbank_recs), genbank_file)
    
    return genbank_recs

#-------------------------------------------------------------------------------
def get_ref_codon(pos_in_gene,seq):
    pos_in_codon = pos_in_gene % 3
    if pos_in_codon==0:
        pos_in_codon=3
    # to make sure we are not in negative numbers with early SNPs
    
    codon_pos = pos_in_gene - pos_in_codon
    if codon_pos < 0:
        
        codon_pos = 0
    ref_codon = str(seq[ codon_pos : (codon_pos + 3)])
    return ref_codon
#-----------------------------------------------------------------------------------------
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
        return "Error non canonical base"
#--------------------------------------------------------------------------------------------
def missing_char(str, pos,n):
    item=list(str)
    item[pos]=n
    str="".join(item)
    return str
#--------------------------------------------------------------------------------------------
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

#--------------------------------End Functions------------------------------------------------------------
options.mode.chained_assignment = None

#-------------------------------Parse arguments-------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', help="combined snp table", default='output.txt')
parser.add_argument('-g', '--genbank', help="genbank file of reference genome")
parser.add_argument('-t', '--translation_table', help="provide translation table number", default=11)
parser.add_argument('-q', '--query', nargs='*', help="list of query files separated by spaces", metavar="qfile.dat")
parser.add_argument('-n', '--name',nargs='*', help="list of galaxy names for galaxy only", default='nfile.txt')

args = parser.parse_args()
output_file = args.output
genbank=args.genbank
t_t=int(args.translation_table)
query=args.query
names=args.name
#LOG_FILENAME='read_snp.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.WARNING)

# -------------------------------concatenates vcf files-------------------------------
data = [read_csv(f, sep ='\t', header=None, dtype=object) for f in query]

#find the genome with the longest list of SNPs
maxim=list()
for d in data:
    maxim.append(d.index.size)

pos=max(xrange(len(maxim)), key=maxim.__getitem__)

#start with the longest list for filling in information
df=data[pos]
for i,x in enumerate(df[4]):
    if x=='.':
        df[4][i]=df[3][i]
    else:
        pass

name_list=[]
for n in names:
    name_list.append(n.split('.t')[0])


#files1=[]
#for f in query:
#    files1.append(f.split('.t')[0])

df.rename(columns={0:'molecule',1:'refpos',3:'refbase',4:('qbase:'+name_list[pos])},inplace=True)

df.drop_duplicates(inplace=True)
df.set_index(['molecule','refpos'], inplace=True)

for v,di in enumerate(data):
    if v==pos:
        pass
    else:
        di.rename(columns={0:'molecule',1:'refpos'}, inplace=True)
        di.drop_duplicates(inplace=True)
        di.set_index(['molecule','refpos'],inplace=True)
        df=concat([df,di.iloc[:,2]],axis=1,join_axes=[df.index])
        for i,x in enumerate(df[4]):
            if x=='.':
                df[4][i]=df['refbase'][i]
            else:
                pass
        df.rename(columns={4:('qbase:'+name_list[v])},inplace=True)

#pdb.set_trace()
df.fillna('No Hit',inplace=True)
df.reset_index(inplace=True)

for i,f in enumerate(df['molecule']):
    if '|' in f:
        df['molecule'][i]=f.split('|')[-1]


#---------------------read in genbank file features-----------------------------------------------
gblist=parse_genbank_file_list(genbank)


#---------------make dataframe with additional columns to fill out-------------------------
filler=repeat('intergenic',df.index.size).tolist()
sec= DataFrame({'gene_name':filler,'gene_start':filler,'gene_end':filler,'gene_length':filler,'snps_per_gene':filler,'pos_in_gene':filler,'ref_codon':filler,'ref_aa':filler,'query_codon':filler,'query_aa':filler,'transition/transversion':filler,'snps/gene_length':filler,'dn/ds':filler,'product':filler,'strand':filler})
sec=sec.reindex_axis(['gene_name','gene_start','gene_end','gene_length','snps_per_gene','pos_in_gene','ref_codon','ref_aa','query_codon','query_aa','transition/transversion','snps/gene_length','dn/ds','product','strand'],axis=1)

df2=concat([df,sec],axis=1)
#pdb.set_trace()


# exact position start in genbank gets normalized but not end. The imported SNPs are not normalized so for calculating the position in the gene add 1 from the start

#--------------------------- get snps in genes & information------------------------------
# gets gene and sequence of SNP hits from genbank file
gene_hits=[]
seqen=[]
for v,f in enumerate(df['molecule']):
    for gb in gblist:
        if gb.name==f:
            for feature in gb.features:
                location=feature.location
                if( feature.type != "source" and int(df['refpos'][v]) >= (location.start+1) and int(df['refpos'][v]) <= location.end ):
                    # here we append the whole gene with the hit
                    gene_hits.append((feature, gb.name, df['refpos'][v]))
                    #append gene sequence
                    if feature.type=='gene':
                        seqen.append((location.extract( gb.seq),feature.qualifiers['locus_tag'][0]))
#pdb.set_trace()

#------------------------------- gene_length, name, start,stop, pos_in_gene, product--------------------------
for i,item in enumerate(df2.index):
    for v,gh in enumerate(gene_hits):
        
        if df2['molecule'][i]==gh[1] and df2['refpos'][i]==gh[2]:
            if gh[0].type=='gene':
                df2['gene_length'][i]=len((gh[0].location))
                df2['gene_name'][i]=gh[0].qualifiers['locus_tag'][0]
                df2['product'][i]='No product'
                if gh[0].strand==1:
                    df2['gene_start'][i]=gh[0].location.start.position
                    df2['gene_end'][i]=gh[0].location.end.position
                    df2['pos_in_gene'][i]=int(gh[2])-(gh[0].location.start.position+1)
                    df2['strand'][i]=1
                else:
                    df2['gene_start'][i]=gh[0].location.end.position
                    df2['gene_end'][i]=gh[0].location.start.position
                    df2['pos_in_gene'][i]=gh[0].location.end.position-int(gh[2])
                    df2['strand'][i]=-1
            else:
                df2['product'][i]=''.join(gh[0].qualifiers['product'])
#pdb.set_trace()

# --------------------reference codon & aa------------------------------
for i,item in enumerate(df2.index):
    for si in seqen:
        if si[1]==df2['gene_name'][i]:
            df2['ref_codon'][i]=get_ref_codon(df2['pos_in_gene'][i],si[0])
            df2['ref_aa'][i]=str(Seq(df2['ref_codon'][i], generic_dna).translate(table=t_t))

#pdb.set_trace()
#---------------get query codon & aa-------------------------------------------
pos1=list()
for i,item in enumerate(df2['pos_in_gene']):
        if item!='intergenic':
            pos_in_codon=item % 3
            pos1.append(pos_in_codon)
        else:
            pos1.append(item)
#pdb.set_trace()
# get query base information

df2.rename(columns={2:'syn?'}, inplace=True)
count_qbase=list(df2.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in v:
        qindexes.append(i)
df3=df2.iloc[:,qindexes]

#pdb.set_trace()
#-------------------------------get allele for each position---------------------------------
snp_nb=list()
for i,item in enumerate(df3.index):
    snp_u=[n for n in unique(df3.iloc[i,:])[:] if n!=df2.refbase[i]]
    if len(snp_u)>0:
        snp_nb.append(snp_u)
    else:
        snp_nb.append([n for n in unique(df3.iloc[i,:])[:]])

for i,item in enumerate(snp_nb):
    if len(item)<1:
        snp_nb[i]=[n for n in unique(df3.iloc[1,:])[:]]
for i,item in enumerate(snp_nb):
    if snp_nb[i][0]==df2.refbase[i]:
        df2['syn?'][i]='No SNP'

#pdb.set_trace()
#-------------------------------transition/transversion -------------------------------
for g,item in enumerate(df2.index):
    tt=list()
    if df2['syn?'][g]=='No SNP' and df2['gene_name'][g]!='intergenic':
        tt.append('No SNP')
    elif df2['gene_name'][g]=='intergenic':
        tt.append('intergenic')
    elif len(snp_nb[g])==1 and df2['gene_name'][g]!='intergenic' and df2['syn?'][g]!='No SNP':
        tt.append(get_trans_state(df2['refbase'][g],snp_nb[g][0]))
    else:
        for n in snp_nb[g]:
            if n!='No Hit':
                tt.append(get_trans_state(df2['refbase'][g],n))
    df2['transition/transversion'][g]=str('/'.join(tt))

#pdb.set_trace()
#-------------------------------get snps per gene-------------------------------
#only genes with SNPs and not identical locations
snps_gene=df2[df2['syn?']!='No SNP'].groupby('gene_name').size().reset_index()
snps_gene2=list()
for i,v in enumerate(snps_gene.index):
    col=[n for n in snps_gene.iloc[i,:][:]]
    snps_gene2.append(col)


for i,v in enumerate(snps_gene2):
    for l,s in enumerate(df2['gene_name']):
        if v[0]==s and s!='intergenic':
            df2['snps_per_gene'][l]=v[1]
        elif s=='intergenic' and df2['syn?'][l]!='No SNP':
            df2['snps_per_gene'][l]='intergenic'
        elif df2['syn?'][l]=='No SNP'and s!='intergenic':
            df2['snps_per_gene'][l]=0



#pdb.set_trace()

#-------------------------------snps per gene length-------------------------------
for x,i in enumerate(df2.snps_per_gene):
    if i==0:
        df2['snps/gene_length'][x]=0
    elif i != 'intergenic':
        df2['snps/gene_length'][x]=float(i)/float(df2.gene_length[x])


# -------------------------------query codon & aa-------------------------------
#pdb.set_trace()

for i,item in enumerate(df2.ref_codon):
    if item != 'intergenic':
        if len(snp_nb[i])==1 and df2['strand'][i]==1:
            item1=missing_char(item,pos1[i],snp_nb[i][0])
            df2['query_codon'][i]=item1
            df2['query_aa'][i]=str(Seq(df2['query_codon'][i], generic_dna).translate(table=t_t))
        elif len(snp_nb[i])==1 and df2['strand'][i]!=1:
            g=invert_nucl(snp_nb[i][0])
            item1=missing_char(item,pos1[i],g)
            df2['query_codon'][i]=item1
            df2['query_aa'][i]=str(Seq(df2['query_codon'][i], generic_dna).translate(table=t_t))
        else:
            mult=list()
            mult_aa=list()
            
            for n in snp_nb[i]:
                if df2['strand'][i]==1 and n!='No Hit':
                    item1=missing_char(item,pos1[i],n)
                    aa1=str(Seq(item1, generic_dna).translate(table=t_t))
                    mult.append(item1)
                    mult_aa.append(aa1)
                elif n!='No Hit':
                    g=invert_nucl(n)
                    item1=missing_char(item,pos1[i],g)
                    aa1=str(Seq(item1, generic_dna).translate(table=t_t))
                    mult.append(item1)
                    mult_aa.append(aa1)
                df2['query_codon'][i]=str('/'.join(mult))
                df2['query_aa'][i]=str('/'.join(mult_aa))

for i,item in enumerate(df2.ref_codon):
    if df2['syn?'][i]=='No SNP':
        df2['query_codon'][i]=item
        df2['query_aa'][i]=df2.ref_aa[i]

# -------------------------------synonymous nonsynonymous-------------------------------
for i,item in enumerate(df2.query_aa):
    if '/' in item:
        mult=list()
        for n in item.split('/'):
            if n==df2.ref_aa[i]:
                mult.append('SYN')
            else:
                mult.append('NSYN')
        df2['syn?'][i]='/'.join(mult)
    else:
        if item=='intergenic':
            df2['syn?'][i]='intergenic'
        elif item==df2.ref_aa[i] and  df2['syn?'][i]!='No SNP':
            df2['syn?'][i]='SYN'
        elif df2['syn?'][i]!='No SNP':
            df2['syn?'][i]='NSYN'

#pdb.set_trace()

#-------------------------------dn/ds-------------------------------
dn=df2.groupby('gene_name')['syn?']
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
df2['dn/ds'].update(Series(flat_dn))
#pdb.set_trace()
#-------------------------------write file-------------------------------
df2.sort(columns=['molecule','refpos'],inplace=True)
#pdb.set_trace()
with open(output_file,'w') as output2:
    df2.to_csv(output2, sep='\t', index=False)


