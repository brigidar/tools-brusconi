#!/usr/bin/env python
#########################################################################################
# Adapted from gbk_to_gff

# Name        :select_product.py
# Version     : 0.2
# Description :Finds all CDS that match the key words provided in the comma separated argument -k in the product description and adds a given number of bases on each side of the transcript for the exclusion set
# Date        :March 10th, 2016
# Author      : Brigida Rusconi
#
# #########################################################################################





import os, re, sys, argparse
import collections
#import pdb
from Bio import SeqIO

def open_file(fname):
    """
        Open the file (supports .gz .bz2) and returns the handler
        
        @args fname: input file name for reading
        @type fname: str
        """
    
    try:
        if os.path.splitext(fname)[1] == ".gz":
            FH = gzip.open(fname, 'rb')
        elif os.path.splitext(fname)[1] == ".bz2":
            FH = bz2.BZ2File(fname, 'rb')
        else:
            FH = open(fname, 'rU')
    except Exception as error:
        sys.exit(error)
    
    return FH


def gbk_parse(fname,keys_1,n,outfile=""):
    """
    Extract genome annotation recods from genbank format 

    @args fname: gbk file name 
    @type fname: str
    """
    fhand = open_file(gbkfname)
    lines = []

    for record in SeqIO.parse(fhand, "genbank"):
        gene_tags = dict()
 
        mol_type, chr_id = None, None 

        for rec in record.features:


            product = "NA"
            
            if rec.qualifiers.has_key('product'):
            
                product = rec.qualifiers['product'][0]


            if rec.type == 'source':
                try:
                    mol_type = rec.qualifiers['mol_type'][0]
                except:
                    mol_type = '.'
                    pass 
                try:
                    chr_id = rec.qualifiers['chromosome'][0]
                except:
                    chr_id = record.name 
                continue 
            

            start = int(rec.location.start)+1-n
            stop = int(rec.location.end)+n

            if start < 0:
                start = 1
                    
            for k in keys_1.split(","):
    
                if str(k) in str(product):
                    lines.append('\t'.join([str(chr_id),str(start), str(stop), product])+"\n")

    
    fhand.close()

    return lines

if __name__=='__main__':

    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('-k', '--keys_1', type=str, help="list of products to key on in a comma delimited list")
        parser.add_argument('-g', '--gbk_file', help="genbank file of reference")
        parser.add_argument('-n', '--number', type=int, help="flanking region to include")
        parser.add_argument('-o', '--output', help="output file")

        args = parser.parse_args()
        gbkfname = args.gbk_file
        n=args.number
        keys_1=args.keys_1
        output=args.output
    
    except Exception,e:
        print "Error: %s" % e
        sys.exit(-1)

    ## extract gbk records

    print "Filtering on keys: %s" % keys_1
    lines = gbk_parse(gbkfname,keys_1,n)

    with open(output, 'w') as of:
        for l in lines:
            of.write(l)




