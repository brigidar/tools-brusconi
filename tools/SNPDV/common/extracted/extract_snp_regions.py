#!/usr/bin/env python

import argparse, os, shutil, sys, tempfile
import pdb

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

"""

Author: Mando Rodriguez, with consulting from Brigida Rusconi

extract snp regions: Extracts a SNP position from the source of a given
genbank file flanked by nucleotides on each side of a given size. Output
is collected into a multifasta file. 

Notes:
    
for intergenics
http://biopython.org/wiki/Intergenic_regions

"""



def parse_snp_positions(snp_file):

    lines = []
    num_snps = 0
    with open(snp_file, "rU") as input_handle:

        for line in input_handle:

            if not line.isspace():

                lines.append( line.split() )
                num_snps += 1

    print "Parsed %d snp positions" % num_snps

    return lines

#-------------------------------------------------------------------------------
# This does the same as above but returns a bucketed python dictionary.
#
#
#-------------------------------------------------------------------------------
def parse_snp_positions_dict(snp_file):

    snps = {}
    num_snps = 0
    with open(snp_file, "rU") as input_handle:

        for line in input_handle:

            if not line.isspace():

                snp_txt = line.split()

                if len(snp_txt) != 2:
                    print "Error parsing entry \"%s\", not in form 'locus position'" % line
                    next

                else:

                    try:
                        snp_pos = int(snp_txt[1])
                    except ValueError:
                        print "Can't parse SNP position %s %s, position isn't a valid value" % (snp_txt[0], snp_txt[1])
                        continue

                    
                    if snps.get(snp_txt[0]) is None:

                        snps[snp_txt[0]] = [snp_pos]

                    else:

                        if not snp_pos in snps[snp_txt[0]]:
                            
                            snps[snp_txt[0]].append(snp_pos)

                        else:

                            print "SNP position %s is already in %s, skipping" % (snp_txt[1], snp_txt[0])
                    
                num_snps += 1

    for key, value in snps.items():

        snps[key].sort()

    print "Parsed %d snp positions" % num_snps

    return snps


#-------------------------------------------------------------------------------

def get_gbank_hits(genbank_recs,molecule_id):
    """
    This will break up the molecule id into separate parts from the pipe, and also
    look for an entry with the '.' in it. If it finds one then it will create a list
    of all parts and use that to determine if there is a hit from the genbank file
    locus tag.
    """

    molecule_id_parts = molecule_id.split('|')

    # probably not needed but I'm going to match on the molecule id as well
    # just in case.
    molecule_id_parts.append(molecule_id)

    for mip in molecule_id_parts:

        if '.' in mip:
            try:
                molecule_id_parts.append( mip.split('.')[0] )
            except Exception:
                print "Error getting '.' prefix from locus tag %s in molecule id %s" % (mip,molcule_id) 

    gbhits = [g for g in genbank_recs if g.name in molecule_id_parts or g.id in molecule_id_parts]

    return gbhits
    
#-------------------------------------------------------------------------------

def __main__():

    #Parse Command Line
    parser = argparse.ArgumentParser()


    parser.add_argument('-f','--flanking-bases',  type=int, default=20, help='Number of bases on *each side* of the SNP positions in the reference genome. Default is 20.' )
    parser.add_argument('-s', '--snp-positions', required=True, help="snp panel file")
    parser.add_argument('-g', '--ref-genbank', required=True, help="genbank annotation file")
    parser.add_argument('-o', '--out', default="extracted_snp_region.txt", help="output file")
    parser.add_argument('-v', '--version', action='version', version='1.3') #version='%(prog)s 1.3')
    

    args = parser.parse_args()

    flanking_bases = int(args.flanking_bases)
    snp_positions_file = args.snp_positions
    ref_genbank = args.ref_genbank
    output_file = args.out

    # First we open the genbank file and load all of the records
    genbank_recs = []

    with open(ref_genbank, "rU") as  genbank_input_handle:

        for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

            genbank_recs.append(gbrec)

    print "Read in %i genbank records" % (len(genbank_recs))


    
    # here we load the snp positions from the panel file into a dict
    snp_dict = parse_snp_positions_dict(snp_positions_file)

    extracted_snps = []

    required_length = flanking_bases + 1

    for molecule_id in snp_dict.keys():

        #
        # Fetch all of the genbank records that match the molecule name
        #
        #gbhits = [g for g in genbank_recs if g.name == molecule_id or g.id == molecule_id]

        gbhits = get_gbank_hits(genbank_recs,molecule_id)

        if len(gbhits) == 0:
            print "Molecule '%s' is not found in genbank file" % molecule_id
            continue

        
        # we go through all of our genbank hits

        for gb in gbhits:

            print "Extracting %d SNP positions for molecule '%s' from genbank '%s'" % (len(snp_dict[molecule_id]), molecule_id, gb.description)

            # we go through all of our positions for each genbank molecule hit.
            for pos in snp_dict[molecule_id]:

                try:
                    
                    # We need to normalize the position "starts" it starts from one
                    # in the genbank but in python list terms it starts at zero.
                    p = pos-1
                        
                    start = p-(flanking_bases)
                    stop = p+(flanking_bases)+1 # add one to this so the last nuc is inclusive

                    seq_end = len(gb.seq)


                    status = 0
                    
                    if stop > seq_end:
                        
                        stop = seq_end

                        status = 1


                    elif start < 0:

                        start = 0

                        status = 2
                        

                    this_seq = gb.seq[start:stop]
                

                    # if it's not of the required length, then we just bail out with an exception and keep going
                    # new required length should at least include the snp. if len(this_seq) != required_length:                        
                    
                    if not len(this_seq) > required_length:                        
                        raise Exception("Not of the required length, is %d should be at least %d" % (len(this_seq), required_length))


                    if status == 1:
                        print "!!!!! Extracted SNP %s is near the end of a sequence and has length %d\n" % ("%s_%s_SUBSEQ" % (molecule_id,pos), len(this_seq))

                    elif status == 2:
                        print "!!!!! Extracted SNP %s is near the beginning of a sequence and has length %d\n" % ("%s_%s_SUBSEQ" % (molecule_id,pos), len(this_seq))

                        

                    
    
                    
                    record = SeqRecord(this_seq,id="%s_%s_SUBSEQ" % (molecule_id,pos),
                                       name="%s_%s_SUBSEQ" % (gb.name,p),
                                       description="Extracted subsequence from %s at position %d" % (molecule_id, pos))

                    extracted_snps.append(record)

                    status = 0
                    
                except Exception, e:

                    print "!Error retrieving SNP region for '%s' at position %d: %s" % (molecule_id, pos, e)



        
    # outside of the loop, now we write them all to a file.

    with open(output_file, "w") as output_handle:
        
        if len(extracted_snps) > 0:
        
            count = SeqIO.write(extracted_snps, output_handle, "fasta")
            print "Wrote out %s extracted snps to '%s'" % (count,output_file)

        else:
        
            output_handle.write("No hits")
            print "No hits for snps in the genbank file\n"

            
if __name__=="__main__": __main__()
