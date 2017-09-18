#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	5_sorting_genotyper2.py								#
# Version     : 1.6									#
# Project     : SNP table downstream analysis						#
# Description : Script to summarize the genotyper output and classify according to groups etc. for tables that use * for stop codon. Transitions and transversions are counted even if in multiallelic state (can be more than total amount of genes, since each allele represents a change.	#
# Author      : Brigida Rusconi								#
# Date        : August 3rd, 2017							#
#											#
#########################################################################################

import argparse, os, sys, csv,pdb
from pandas import *

#-------------------argparse arguments-----------------------------------------------------------------------------------------

#output and input file name
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="genotype group summary")
parser.add_argument('-s', '--genotyper_table', help="genotyper table to summarize", default="snp_genotype_out.txt")

parser.add_argument('-r', '--refgenome', help="refgenome name")

args = parser.parse_args()
output_file = args.output
input_file = args.genotyper_table
ref_genome= args.refgenome

#----------------file reading--------------------------------------------------------------------------------------------


#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0, dtype=str)


#-------------------------PI/NI counting-----------------------------------------------------------------------------------

##get sums of total patterns, SYN, NSYN, Intergenic according to PI, NI and multiallelic state

# total snp
total=[df.index.size]

#total of genes with snp
genes=[(df.groupby('gene_name').size().index.size)-1]

#total PI/NI
nl_t=[df[df['Informative']=='NI'].index.size]
pl_t=[total[0]-nl_t[0]]
# get total of PI NI of genes
gene_inf=df.groupby(['gene_name','Informative']).size().reset_index()
nl=[gene_inf[(gene_inf['Informative']=='NI')].index.size]
pl=[gene_inf[(gene_inf['Informative']!='NI')].index.size]

#--------------------------------------------stop count----------------------------------------------------------------

#total stop gains
stop_t=[df[(df['syn?'].str.contains('NSYN')) & (df['query_aa'].str.contains('\*'))].index.size]

#genes with gain of stop
st=df.groupby(['gene_name','syn?','query_aa']).size().reset_index()
stop_g=[st[(st['syn?'].str.contains('NSYN')) & (st['query_aa'].str.contains('\*'))].index.size]

#total stop losses
stop2_t=[df[(df['syn?'].str.contains('NSYN')) & (df['ref_aa'].str.contains('\*'))].index.size]
                                   
#genes with loss of stop
st2=df.groupby(['gene_name','syn?','ref_aa']).size().reset_index()
stop2_g=[st[(st['syn?'].str.contains('NSYN')) & (st2['ref_aa'].str.contains('\*'))].index.size]

#--------------------------------------hypothetical count----------------------------------------------------------------------
#hypothetical SNPs
hyp_t=[df[df['product'].str.contains('hypothetical')].index.size]
                                   
#hypothetical genes
hyp=df.groupby(['gene_name','product']).size().reset_index(level=1)
hypo=[hyp[hyp['product'].str.contains('hypothetical')].index.size]

#------------------------------------------multiallelic count------------------------------------------------------------------

#total amount of multiallelic positions
multi_t=[(df['syn?'].str.contains('/').sum())]
                                   
#multallelic genes
mult=df.groupby(['gene_name','query_codon']).size().reset_index()
multi=[mult[mult['query_codon'].str.contains('/')].index.size]

#-------------------------------------------transitions/transversions count-----------------------------------------------------------------
#If more than 2 alleles will count total amount of transitions and transversions, which will be higher if there are multiallelic positions.

#total transitions
trans_t=[df[df['transition/transversion']=='transition'].index.size]
#total transversions
trvs_t=[df[df['transition/transversion']=='transversion'].index.size]

for i,v in enumerate(df['transition/transversion']):
    ts=trans_t[0]
    tv=trvs_t[0]
    if '/' in v:
        m=v.split('/')
        ts+=m.count('transition')
        tv+=m.count('transversion')
    trans_t[0]=ts
    trvs_t[0]=tv


gene_pro=df.groupby(['gene_name','transition/transversion']).size().reset_index()
#Gene transitions
trans_g=[gene_pro[gene_pro['transition/transversion'].str.contains('transition')].index.size]

#Gene transversions
trvs_g=[gene_pro[gene_pro['transition/transversion'].str.contains('transversion')].index.size]

#------------------------------------------total columns------------------------------------------------------------------
#first row total values
tot=total+nl_t+pl_t+stop_t+stop2_t+hyp_t+trans_t+trvs_t+multi_t
 #genes with column
t_genes=genes+nl+pl+stop_g+stop2_g+hypo+trans_g+trvs_g+multi

#---------------------------------------------subsetting by snp type---------------------------------------------------------------
#counting number of syn/nsyn
data=DataFrame(df.groupby('syn?').size())


#-------------------------------------------------PI/NI by subtype-----------------------------------------------------------
#distribution of informative non-informative
pi_dis= DataFrame(df.groupby(['syn?', 'Informative']).size()).reset_index(level=1)
inf=DataFrame(index=pi_dis.index)
for i,group in pi_dis.groupby('Informative'):
    inf=concat([inf,group[0]], axis=1)
data=concat([data,inf], axis=1)

#------------------------------------------------stop by subtype------------------------------------------------------------
# distribution of stops query, remove synonymous changes

si_dis= DataFrame(df.groupby(['syn?','query_aa']).size()).reset_index()

tes=si_dis[si_dis['query_aa'].str.contains('\*')]
si_dis2=tes.reset_index().groupby('syn?').sum()
for i,v in enumerate(si_dis2.index):
    if 'NSYN' not in v:
        si_dis2[0][i]=0

data=concat([data,si_dis2.iloc[:,1]], axis=1, join_axes=[data.index])

# distribution of stops reference, remove synonymous changes
si_dis3= DataFrame(df.groupby(['syn?','ref_aa']).size()).reset_index()
tes3=si_dis3[si_dis3['ref_aa'].str.contains('\*')]
si_dis3=tes3.reset_index().groupby('syn?').sum()
for i,v in enumerate(si_dis3.index):
    if 'NSYN' not in v:
        si_dis3[0][i]=0

data=concat([data,si_dis3.iloc[:,1]], axis=1, join_axes=[data.index])

#-------------------------------------hypothetical by subtype-----------------------------------------------------------------------
# distribution of hypothetical
hy_dis= DataFrame(df.groupby(['product', 'syn?']).size()).reset_index(level=0)
hy2=hy_dis[hy_dis['product'].str.contains('hypothetical')]
hy3=hy2.reset_index().groupby('syn?').sum()
data=concat([data,hy3], axis=1, join_axes=[data.index])

#--------------------------------------transitions/transversions by subtype----------------------------------------------------------------------
#distribution of transition transversion
ti_dis=DataFrame(df.groupby(['syn?', 'transition/transversion']).size()).reset_index()
ti=DataFrame(index=ti_dis.index)
tl=[]
r=[]
rl=[]
#get one value for transition and one for transversion
for i,group in ti_dis.groupby('syn?'):
    tl.append([n for n in group['transition/transversion'].values])
    rl.append([n for n in group.iloc[:,2].values])#[n for n in
trans=[0 for i in range(0,len(rl))]
trasv=[0 for i in range(0,len(rl))]
for i,m in enumerate(tl):
    for t,c in enumerate(m):
        p=c.split('/')
        for d in p:
            if d=='transition':
                trans[i]+= rl[i][t]
            elif d=='transversion':
                trasv[i]+= rl[i][t]
data.insert(data.columns.size, 'transition', trans)
data.insert(data.columns.size, 'transversion', trasv)

#---------------------------------------------multiallelic by subtype---------------------------------------------------------------
# distribution of multiallelic
mu_dis= DataFrame(df.groupby(['syn?']).size()).reset_index()
for i,v in enumerate(mu_dis['syn?']):
    if '/' not in v:
        mu_dis.iloc[i,1]=0
mu_dis=mu_dis.rename(columns = {0:'mul'}).set_index(['syn?'])

#-----------------------------------------------concatenating table-------------------------------------------------------------
data=concat([data,mu_dis['mul']], axis=1, join_axes=[data.index])
data=data.T
data.insert(0,'total',tot)
data.insert(1,'Genes_with',t_genes)
col=['SNPs']+['NI']+['PI']+['stop_gain']+['stop_loss']+['hypothetical proteins']+['transition']+['transversion']+['multiallelic']
data.insert(0, 'header',col)
data = data.fillna('0').set_index('header').astype('float64')

#-------------------------------------------------summary genic-----------------------------------------------------------
#add a summary genic column by selecting all columns that have SNPs in genes
hed=data.columns[data.columns.str.contains('SYN')]
gen=data.loc[:,hed]
genic=[]
for i in range(0,gen.index.size):
    g=gen.iloc[i,:].sum()
    genic.append(g)
data.insert(data.columns.size,'genic',genic)

#-------------------------------------------------save file-----------------------------------------------------------
#save summary table
with open(output_file ,'w') as output:
    data.to_csv(output, sep='\t')

#---------------------------------------------transformation to percentage---------------------------------------------------------------
#calculate percentage
perc=data.copy(deep=True)
perc2=DataFrame(index=[data.index])
#calculate percentage for each value
for i in range(0,data.columns.size):
    mol=[n for n in (data.iloc[:,i] / data.iloc[0,i]  * 100)]
    mol1=[n for n in (data.iloc[-3:-1,i] / (data.iloc[0,i] + data.iloc[-1,i])  * 100)]
    mol2=[data.iloc[0,i]/data.loc[data.index[0],'total']*100]
    mol1.append(mol[-1])
    mol3=mol2+mol[1:-3]+mol1
    perc2.insert(i,('% '+ str(data.columns[i])),mol3)

#----------------------------------------------save percentage file--------------------------------------------------------------

#save summary only percentage
with open("percentage_summary.txt" ,'w') as output:
    perc2.to_csv(output, sep='\t')

#----------------------------------------------summary file--------------------------------------------------------------
#summary of all groups and positions and genome names #http://wesmckinney.com/blog/filtering-out-duplicate-dataframe-rows/
grouped = df.groupby('Group')
index = [gp_keys[0] for gp_keys in grouped.groups.values()]
unique = df.reindex(index).iloc[:,3:12].set_index('Group')
#groups the amount of syn nsyn intergenic for each one of the groups
ls1=df.groupby(['Group', 'syn?']).size()
#keeps only the group as an index
ls2=ls1.reset_index(level=1)
#change the name of the counting column
ls2.rename(columns={0:'count_snp'}, inplace=True)
#groups the amount of PI/NI for each one of the groups after grouping by SYN/NSYN
lsp=df.groupby(['Group', 'syn?','Informative']).size()
#keeps only the group as an index
lsp2=lsp.reset_index(level=[1,2])
#change the name of the counting column
lsp2.rename(columns={0:'count_PI'}, inplace=True)

#------------------------------------------------------------------------------------------------------------
#take the refpos that match that group and put them with the molecule names
blist=[]
grouped=df.groupby(['Group', 'syn?'])
for name, group in grouped:
    #restart the list for each subgrouping according to the type of mutation
    mol_list=[]
    for a,b in group.groupby(['molecule']):
        rf=[str(n) for n in b['refpos']]
        mol_list.append(a+ ':'+'/'.join(rf))

    blist.append(', '.join(mol_list))
#concatenates based on the longer index with the repeated groups with the join_axes command
lst=concat([ls2, unique],axis=1, join_axes=[ls2.index])


#------------------------------------------------------------------------------------------------------------
#add column for group that is specific to the reference genome. if the group has only one kind of mutation syn/intergenic, etc then it will call a string instead of a series so we need to check first the class of the loc calling. takes argument with ref genome name.
ni_ref=[]
for i in lst.index:
    m=[]
    pat = ""
    
    if isinstance(lst.loc[i, 'Pattern'], str):
    
        pat = lst.loc[i, 'Pattern']
    
    else:
    
        pat = lst.loc[i, 'Pattern'][0]
    
    
    for n in pat:
        m.append(int(n))

    if all(x==1 for x in m[1:]):
        ni_ref.append(ref_genome)
    else:
        ni_ref.append('--')

#------------------------------------------------------------------------------------------------------------
#remove old PI column
lst.insert(5, 'refpos split up',blist)
lst.insert(lst.columns.size,'ref_genome', ni_ref)

#----------------------------------------------save summary file--------------------------------------------------------------
sum_g=concat([lst.iloc[:,0:2],lst.iloc[:,4:]], axis=1)

#save summary for group
with open("summary_groups.txt",'w') as output:
    sum_g.to_csv(output, sep='\t')












