#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	vcf_snp_output.py								#
# Version     : 0.4								#
# Project     : snp verify reads						#
# Description : Script to populate SNPs identified from reads with location information		#
# Author      : Brigida Rusconi								#
# Date        : September 13th, 2017						#
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
import itertools
from itertools import chain
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

#-------------------------------------------------------------------------------

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
        return "--"
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

#--------------------------------------------------------------------------------------------
def flattern(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flattern(i))
        else: rt.append(i)
    return rt

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

name_list=[]
for n in names:
    name_list.append(n.split('.t')[0])

#start reading in vcf as tables
#https://pandas.pydata.org/pandas-docs/stable/merging.html#concatenating-using-append

df1=data[0]
df1.rename(columns={0:'molecule',1:'refpos',2:'s',3:'refbase',4:('qbase:'+name_list[0])},inplace=True)
df1=df1.iloc[:,0:5]

for i,d in enumerate(data):
    if i > 0:
        pl=data[i].rename(columns={0:'molecule',1:'refpos',2:'s',3:'refbase',4:('qbase:'+name_list[i])})
        df1=df1.merge(pl.iloc[:,0:5], on=['molecule','refpos','refbase','s'],how='outer',sort=False)

df1.drop('s',axis=1,inplace=True)
print "Read vcf files"

for i,f in enumerate(df1['molecule']):
    try:
        if '|' in df1['molecule'][i]:
            df1['molecule'][i]=f.split('|')[-1]
    except KeyError:
        continue


#qindexes w/o molecule and refpos as indexes!!!
count_qbase=list(df1.columns.values)
qindexes=[]
for i, v in enumerate(count_qbase):
    if 'qbase:' in str(v):
        qindexes.append(i)
df1.fillna('No Hit', inplace=True)
df1.mask(df1['refbase'].str.len()>1,inplace=True)
df1.dropna(inplace=True)
#check if any complex event are still present, should be obsolete with new freebayes settings
comple=df1.iloc[:,qindexes]
colm=list(comple.columns.values)
for i in colm:
    comple.mask(comple[i].astype(str).str.contains(','),inplace=True)
clean_df=comple.dropna()
df1=df1[df1.index.isin(clean_df.index)]




#count snps per genome
counting=[]
for i in qindexes:
    b=df1.iloc[:,i].str.contains('\.').sum()
    c=df1.iloc[:,i].str.contains('No Hit').sum()
    counting.append((df1.index.size-(b+c)))

snp_genome=DataFrame({'genomes':name_list, 'snps':counting})

with open('snps_per_genome.txt','w') as output3:
    snp_genome.to_csv(output3,index=False,sep='\t')
#replaces all identical hits with actual nucleotide
for v in df1.index:
#df1.loc[v]=[df1.refbase[v] for x in df1.loc[v] if x=='.']
    df1.loc[v].replace({'\.':df1.refbase[v]},regex=True,inplace=True)


del(data)

#splits up molecule name to match genbank record

df=df1.drop_duplicates(subset=['molecule','refpos'])
#---------------------read in genbank file features-----------------------------------------------

gblist=parse_genbank_file_list(genbank)



# exact position in genbank gets normalized. The imported SNPs are not normalized (start @ 0 instead of 1) use refpos_norm instead

df1['refpos']=df1['refpos'].astype(int)
df1['refpos_norm']=df1['refpos']-1

pl=df1.groupby('molecule')['refpos_norm']

seqen=dict()
gen2=[]
ml=[]
pos2=[]
st2=dict()
sp2=dict()
str2=dict()
prod=dict()
tag2=[]
len2=dict()

for n,g in pl:
    gen=[]
    mol=[]
    pos=[]
    tag1=[]
    for gb in gblist:

        if gb.name==n:
           
            for feature in gb.features:
                if  feature.type != "source":
                    location=feature.location
                    g=g.astype(int)
                    gl=array(g)
                    #the last number is not included in the actual gene. upper limits are never included they are up to references e.g. [2:3] gives nucleotide in position 2 only
                    pr=gl[(gl>=(location.start.position))*(gl< location.end.position)]
                    
                    
                    if pr.size>0:
                        sf=[]
                        yt=[]
                        tag=[]
                        
                        
                        if feature.type=='gene':
                            #get feature info, gene length, molecule, refpos
                            sf=repeat(feature, pr.size).tolist()
                            yt=repeat(gb.name, pr.size).tolist()
                            len2[feature.qualifiers['locus_tag'][0]]=len(location)
                            #get sequence and gene_name
                            seqen[feature.qualifiers['locus_tag'][0]]=location.extract(gb.seq)
                            tag=repeat(feature.qualifiers['locus_tag'][0],pr.size).tolist()
                    
                            if feature.strand==1:
                                #get start stop and strand
                                st2[feature.qualifiers['locus_tag'][0]]=location.start.position
                                sp2[feature.qualifiers['locus_tag'][0]]=location.end.position
                                str2[feature.qualifiers['locus_tag'][0]]='1'
                            else:
                                sp2[feature.qualifiers['locus_tag'][0]]=location.start.position
                                st2[feature.qualifiers['locus_tag'][0]]=location.end.position
                                str2[feature.qualifiers['locus_tag'][0]]='-1'
                            
                            pos.append(pr.tolist())
                            gen.append(sf)
                            mol.append(yt)
                            tag1.append(tag)
                    
                        elif feature.type=='CDS':
                            # make dictionary for product
                            prod[feature.qualifiers['locus_tag'][0]]=''.join(feature.qualifiers['product'])
    gen2.append(gen)
    ml.append(mol)
    pos2.append(pos)
    tag2.append(tag1)
#----------------------------------------------------------

# flatten list of lists with strings https://stackoverflow.com/questions/17864466/flatten-a-list-of-strings-and-lists-of-strings-and-lists-in-python
ml=flattern(ml)
gen2=flattern(gen2)
pos2=flattern(pos2)
tag2=flattern(tag2)
#make table from lists and dictionaries
info=[ml,pos2,tag2]
gene_hits=[ml,pos2,gen2]
table1=DataFrame(info)
table1=table1.T
table1.columns=['molecule','refpos_norm','gene_name']

#map dictionary to column https://stackoverflow.com/questions/24216425/adding-a-new-pandas-column-with-mapped-value-from-a-dictionary
table1['gene_start']=table1['gene_name'].map(st2)
table1['gene_end']=table1['gene_name'].map(sp2)
table1['gene_length']=table1['gene_name'].map(len2)
table1['strand']=table1['gene_name'].map(str2)
table1['product']=table1['gene_name'].map(prod)
table1['pos_in_gene']=table1['refpos_norm']-table1['gene_start']+1
table1['product'].fillna('No CDS',inplace=True)

#replace inverted gene position
for i in table1.index:
    if table1['pos_in_gene'][i]<=0:
        table1['pos_in_gene'][i]=table1['gene_start'][i]-table1['refpos_norm'][i]


##only genes with SNPs and not identical locations
snps_gene=table1.groupby('gene_name').size().reset_index()
snps_gene2=dict(zip(snps_gene['gene_name'].tolist(),snps_gene[0].tolist()))
table1['snps_per_gene']=table1['gene_name'].map(snps_gene2)
table1['snps/gene_length']=table1['snps_per_gene'].astype(float)/table1['gene_length'].astype(float)
table1['refpos']=table1['refpos_norm']+1
print " Read gene information"


# --------------------reference codon & aa------------------------------
#only get reference codon and amino acids for coding genes

table2=table1[table1['product']!='No CDS']
genes=table1[table1['product']=='No CDS']

pos_in_gene=table2['pos_in_gene'].tolist()
tag3=table2['gene_name'].tolist()
ref_codon=[]
for i in range(0,table2.index.size):
    ref_codon.append(get_ref_codon(pos_in_gene[i],seqen[tag3[i]]))


ref_aa=[]
for i in range(0,table2.index.size):
    ref_aa.append(str(Seq(ref_codon[i], generic_dna).translate(table=t_t)))

table2.insert(table2.columns.size,'ref_codon',ref_codon)
table2.insert(table2.columns.size,'ref_aa',ref_aa)
print "Read ref codon and aa"




#SNPs in multiple genes
fl2=table2.append(genes)
tb=fl2.drop_duplicates(subset=['molecule','refpos'])

duplicates=table1[~table1.index.isin(tb.index)]

print('%s SNPs are located in more than one gene') % duplicates.index.size
#---------------get query codon & aa-------------------------------------------

#get alleles to identify query codon and amino acid
tb.set_index(['molecule','refpos'],inplace=True)
df1.set_index(['molecule','refpos'],inplace=True)
coding=df1[df1.index.isin(tb.index)]
coding.reset_index(inplace=True)

# get query base information
df3=coding.iloc[:,qindexes].join(coding.refbase)
#position in codon
pos1=(tb.pos_in_gene.astype(int) % 3).tolist()
pos1=[ x if x!=0 else 3 for x in pos1 ]
pos1=[(x-1) for x in pos1]
ref_codon=tb.ref_codon.astype(str).tolist()
#get allele for each position
snp_nb, idn =get_snp(df3)




query_codon=[]
for i,v in enumerate(snp_nb):
    if ref_codon[i]=='nan':
        query_codon.append('nan')
    else:
        #positions on -1 strand need to be inverted use invert_nucl function
        if tb['strand'][i]==1:
            ts=list()
            if len(v)==1 and v[0]!='No Hit':
                query_codon.append(missing_char(str(ref_codon[i]),pos1[i],v[0]))
            else:
            #multiallelic position gets codon for each
                for n in v :
                    if n!='No Hit':
                        ts.append(missing_char(str(ref_codon[i]),pos1[i],n[0]))
                query_codon.append('/'.join(ts))
        else:
            ts=list()
            if len(v)==1 and v[0]!='No Hit':
                query_codon.append(missing_char(str(ref_codon[i]),pos1[i],invert_nucl(v[0])))
            else:
                #multiallelic position gets codon for each
                for n in v:
                    if n!='No Hit':
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

tb.insert(tb.columns.size,'query_codon',query_codon)
tb.insert(tb.columns.size,'query_aa',query_aa)

print "Read query codons and aa"
#-------------------------------transition/transversion -------------------------------
df1.reset_index(inplace=True)
#df1.sort_values(['molecule','refpos'],inplace=True)
#
df4=df1.iloc[:,qindexes].join(df1.refbase)


snp_nb2, ident=get_snp(df4)
ts_tv=[]
for i,v in enumerate(snp_nb2):
    ts=[]
    if len(v)==1 and v[0]!= 'No Hit':
        ts_tv.append(get_trans_state(df1.refbase[i],v[0]))
    else:
        for n in v:
            if n!='No Hit':
                     ts.append(get_trans_state(df1.refbase[i],n))
        ts_tv.append('/'.join(ts))
df1.insert(df1.columns.size,'transition/transversion',ts_tv)
print "Read transition/transversion"
df1.drop('refpos_norm',axis=1,inplace=True)
del(table1)
del(table2)
del(df3)
del(df4)

df1.set_index(['molecule','refpos'],inplace=True)

fin=df1.join(tb)
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
#-------------------------------write file-------------------------------
with open(output_file,'w') as output2:
    final.to_csv(output2, sep='\t', index=False)


#pdb.set_trace()
#with open('snp_multiple_genes.txt','w') as output2:
#    duplicates.to_csv(output2, sep='\t', index=False)




