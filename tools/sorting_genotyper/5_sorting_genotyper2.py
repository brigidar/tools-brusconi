#!/usr/bin/env python
#########################################################################################
#											#
# Name	      :	5_sorting_genotyper2.py								#
# Version     : 1.6									#
# Project     : SNP table downstream analysis						#
# Description : Script to summarize the genotyper output and classify according to groups etc. for tables that use * instead of Stop for stop codons		#
# Author      : Brigida Rusconi								#
# Date        : March 10th, 2016							#
#											#
#########################################################################################

import argparse, os, sys, csv
import pdb
from pandas import *

#------------------------------------------------------------------------------------------------------------

#output and input file name
parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="genotype group summary")
parser.add_argument('-s', '--genotyper_table', help="genotyper table to summarize", default="snp_genotype_out.txt")

parser.add_argument('-r', '--refgenome', help="refgenome name")

args = parser.parse_args()
output_file = args.output
input_file = args.genotyper_table
ref_genome= args.refgenome

#------------------------------------------------------------------------------------------------------------


#read in files as dataframe
df =read_csv(input_file, sep='\t', header=0, dtype=str)


#------------------------------------------------------------------------------------------------------------

##get sums of total patterns, SYN, NSYN, Intergenic according to PI, NI and multiallelic state

total=[df.index.size]


# get total of PI NI and add to total snps
pl=[]
nl=[]
for name, group in df.groupby('Informative'):
    nl.append(name)

for i in df.groupby('Informative').size():
    pl.append(i)
#------------------------------------------------------------------------------------------------------------


#amount of genes with gain or kept stop
st=df.groupby(['gene_name','query_aa']).size().reset_index(level=1)
stop=[(st['query_aa'].str.contains('\*').sum())]

#amount of genes with loss or kept stop
st2=df.groupby(['gene_name','ref_aa']).size().reset_index(level=1)
stop2=[(st2['ref_aa'].str.contains('\*').sum())]
#------------------------------------------------------------------------------------------------------------

#amount of hypothetical genes
hyp=df.groupby(['gene_name','product']).size().reset_index(level=1)
hypo=[(hyp['product'].str.contains('hypothetical').sum())]
#------------------------------------------------------------------------------------------------------------
#total amount of genes

count=[]
for i,v in enumerate(df.groupby('gene_name').size()):
    count.append(i)
genes=[count[-1]]

#total amount of multiallelic positions
multi=[(df['query_codon'].str.contains('/').sum())]
#------------------------------------------------------------------------------------------------------------

# counting number of transitions transversions:

ti_n=[]
ti_g=[]
for name, group in df.groupby('transition/transversion'):
    ti_n.append(name)
for i in df.groupby('transition/transversion').size():
    ti_g.append(i)

# make list with transition and transversion only to add to it the others with multiallelic state

ti_n2=[]
ti_g2=[]
for i,v in enumerate(ti_n):
    m=v.split('/')
    if len(m)==1:
        ti_n2.append(v)
        ti_g2.append(ti_g[i])



for i,v in enumerate(ti_n):
    m=v.split('/')
    if len(m)>1:
        for t in m:
            if t==ti_n2[0]:
                ti_g2[0]=ti_g2[0] + ti_g[i]
            elif t==ti_n2[1]:
                ti_g2[1]=ti_g2[1] + ti_g[i]


if len(ti_n2)>2:
    sys.exit("error in transition/transversion column check snp table")


#------------------------------------------------------------------------------------------------------------
#first row total values
tot=total+pl+genes+stop+stop2+hypo+ti_g2+multi
#------------------------------------------------------------------------------------------------------------

#counting number of syn/nsyn
other_syn=DataFrame(df.groupby('syn?').size())


#------------------------------------------------------------------------------------------------------------

#distribution of informative non-informative
pi_dis= DataFrame(df.groupby(['syn?', 'Informative']).size()).reset_index(level=1)
test=DataFrame(index=pi_dis.index)
for i,group in pi_dis.groupby('Informative'):
    test=concat([test,group[0]], axis=1)
test2=concat([other_syn,test], axis=1)

#------------------------------------------------------------------------------------------------------------
#distribution of genes:
gi_dis=df.groupby(['syn?', 'gene_name']).size().reset_index(level=0)
gi_dis2=gi_dis.groupby('syn?').size()
tes2=concat([test2,gi_dis2], axis=1)

#------------------------------------------------------------------------------------------------------------

# distribution of stops query
si_dis= DataFrame(df.groupby(['syn?','query_aa']).size()).reset_index()
tes=si_dis[si_dis['query_aa'].str.contains('\*')]
si_dis2=tes.reset_index().groupby('syn?').sum()

test3=concat([tes2,si_dis2.iloc[:,1]], axis=1, join_axes=[test.index])
# distribution of stops reference
si_dis3= DataFrame(df.groupby(['syn?','ref_aa']).size()).reset_index()
tes3=si_dis3[si_dis3['ref_aa'].str.contains('\*')]
si_dis3=tes3.reset_index().groupby('syn?').sum()

tes3=concat([test3,si_dis3.iloc[:,1]], axis=1, join_axes=[test.index])

#------------------------------------------------------------------------------------------------------------


# distribution of hypothetical
hy_dis= DataFrame(df.groupby(['product', 'syn?']).size()).reset_index(level=0)
hy2=hy_dis[hy_dis['product'].str.contains('hypothetical')]
hy3=hy2.reset_index().groupby('syn?').sum()
test4=concat([tes3,hy3], axis=1, join_axes=[test.index])
#------------------------------------------------------------------------------------------------------------


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



test4.insert(test4.columns.size, 'transition', trans)
test4.insert(test4.columns.size, 'transversion', trasv)




#------------------------------------------------------------------------------------------------------------


# distribution of multiallelic
mu_dis= DataFrame(df.groupby(['query_codon', 'syn?']).size()).reset_index(level=0)
mu2=mu_dis[mu_dis['query_codon'].str.contains('/')]
mu3=mu2.reset_index().groupby('syn?').sum()


#get ditsribution of the multiallelic states in the intergenic regions
mu_int=df.groupby('syn?').get_group('intergenic')
counter=0
for i,n in enumerate(mu_int['snp_total']):
    if len(n)==2:
        counter+=1
    elif len(n)==3:
        counter+=2
mu4=mu3.T
mu4.insert(mu4.columns.size, 'intergenic', counter)
mu5=mu4.T
#------------------------------------------------------------------------------------------------------------
test5=concat([test4,mu5], axis=1, join_axes=[test.index])
test5=test5.T
test5.insert(0,'total',tot)
col=['SNPs']+[n for n in nl]+['genes']+['stop_gain']+['stop_loss']+['hypothetical proteins']+[n for n in ti_n2]+['multiallelic']
test5.insert(0, 'header',col)
test6 = test5.fillna('0').set_index('header').astype('float64')
#remove count of intergenic in gene name
test6.loc['genes','intergenic']=0
#remove syn stop mutations
names=[]
for i in test6.columns.values:
    names.append(i)
if 'SYN' in names:
    test6.loc['stop_gain','total']=test6.loc['stop_gain','total']-test6.loc['stop_gain','SYN']
    test6.loc['stop_loss','total']=test6.loc['stop_loss','total']-test6.loc['stop_loss','SYN']
    test6.loc['stop_loss','SYN']=test6.loc['stop_loss','SYN']-test6.loc['stop_loss','SYN']
    test6.loc['stop_gain','SYN']=test6.loc['stop_gain','SYN']-test6.loc['stop_gain','SYN']



# check back to string

#------------------------------------------------------------------------------------------------------------
#add a summary genic column
genic=[]
for i,v in enumerate(test6.index):
    gen=[n for n in test6.iloc[i,1:-1][:]]
    genic.append(sum(gen))


test6.insert(test6.columns.size,'genic',genic)
#------------------------------------------------------------------------------------------------------------
#pdb.set_trace()
#calculate percentage
perc=test6.copy(deep=True)
perc2=DataFrame(index=[test6.index])
#calculate percentage for each value
for i in range(0,test6.columns.size):
    mol=[n for n in (test6.iloc[:,i] / test6.iloc[0,i]  * 100)]
    # since transition and transversion are separate for the multiallelic states the total has to be double for two state does not account for 3 state yet.
    mol1=[n for n in (test6.iloc[-3:-1,i] / (test6.iloc[0,i] + test6.iloc[-1,i])  * 100)]
    mol2=[test6.iloc[0,i]/test6.loc[test6.index[0],'total']*100]
    mol1.append(mol[-1])
    mol3=mol2+mol[1:-3]+mol1
    perc2.insert(i,('% '+ str(test6.columns[i])),mol3)
#pdb.set_trace()
#------------------------------------------------------------------------------------------------------------

#save summary table
with open(output_file ,'w') as output:
    perc.to_csv(output, sep='\t')
#------------------------------------------------------------------------------------------------------------

#save summary only percentage
with open("percentage_summary.txt" ,'w') as output:
    perc2.to_csv(output, sep='\t')
#------------------------------------------------------------------------------------------------------------

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
#groups the amount of PI?NI for each one of the groups after dividing by SYN/NSYN
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

ni_ref=[]
#------------------------------------------------------------------------------------------------------------


#add column for group that is specific to the reference genome. if the group has only one kind of mutation syn/intergenic, etc then it will call a string instead of a series so we need to check first the class of the loc calling. takes argument with ref genome name.
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
#------------------------------------------------------------------------------------------------------------

sum_g=concat([lst.iloc[:,0:2],lst.iloc[:,4:]], axis=1)
#pdb.set_trace()
#save summary for group
with open("summary_groups.txt",'w') as output:
    sum_g.to_csv(output, sep='\t')












