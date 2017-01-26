



import Bio
from Bio import Entrez
from Bio import SeqIO
import argparse, os, sys
#import pdb

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genomes', help="ncbi id")
parser.add_argument('-e','--emails', help='email address')

args = parser.parse_args()
genomes=args.genomes
emails=args.emails

Entrez.email= emails
records= Entrez.efetch(db='nucleotide',id=genomes,rettype="gbwithparts",retmode="text")
with open("output.gb",'w') as output2:
    output2.write(records.read())

genomes = SeqIO.read("output.gb", "genbank")
#header = ">"+ genomes.id + " " + genomes.features[0].qualifiers['organism'][0]
#sequence = genomes.seq()

with open("output.fasta",'w') as output2:
    output2.write(genomes.format("fasta"))