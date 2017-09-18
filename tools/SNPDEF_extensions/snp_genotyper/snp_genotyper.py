#!/usr/bin/env python

import argparse

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import ExtendedIUPACDNA
import pdb

#*******************************************************************************
# Program: SNP Genotyper
#
# Author: Mando Rodriguez, with consulting from Brigida Rusconi
#
#
# Based on the original code here:
#   https://github.com/Henshina/Genotyper/blob/master/Geno1.py
#*******************************************************************************


#*******************************************************************************
# This function is used to translate specific codons to amino acids based on
# different translation tables. If a table isn't given then it just does a translation
# using the IUPAC ambiguous DNA table in biopython.
#
# Taken from the ncbi genetics codes here:
#   http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2
#
#   http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG11
#
#
# Some info on the biopython alphabets, may need to use a different one for different
# cases.
#
# http://biopython.org/DIST/docs/api/Bio.Alphabet-module.html
#
# http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC-module.html
#
#*******************************************************************************
def translate_codon(codon, table=1):
    if codon == "--":
        return "--"
    elif codon == "indel":
        return "indel"
    return str(Seq(codon, ExtendedIUPACDNA()).translate(table=table))

#************* end translate_codon function *********************





#*** SNP class object ******************************************
# Bulding functionality into this class object
# to save a bunch of work.
#***************************************************************
class SNP:
    
    #--------------------------------------------
    def __init__(self, qindexes, snp_row, query_genomes, indexes, table=1):

        self.snp_data = snp_row.split('\t')

        # molecule is at index 0
        self.molecule = self.snp_data[0]
        
        self.qindexes = qindexes

        ref_pos_index = indexes['refpos']
        ref_base_index = indexes['refbase']
        ref_codon_index = indexes['ref_codon']
        ref_aa_index = indexes['ref_aa']
        gene_name_index = indexes['gene_name']
        query_codon_index = indexes['query_codon']
        query_aa_index = indexes['query_aa']

        # the ref pos
        self.ref_pos = self.snp_data[ref_pos_index]

        # just the ref base
        self.ref_base = self.snp_data[ref_base_index]

        # collect the qbases from the SNP
        self.qbases = [self.snp_data[qi] for qi in self.qindexes]

        self.query_genomes = query_genomes

        # get the ref counts from the SNP
        self.pattern_list = self.get_pattern()

        # Convert our ref counts to a string as an ID for hashing.
        self.pattern = ''.join([str(b) for b in self.pattern_list])

        # the occurance of different numbers in the pattern
        self.zeroes = self.pattern_list.count(0)
        self.ones = self.pattern_list.count(1)
        self.twos = self.pattern_list.count(2)
        self.threes = self.pattern_list.count(3)

        self.info = self.set_informative()


        # lists of the genomes at the positions with the number in the pattern
        self.genomes_w_zero = self.list_to_string(self.get_genomes(0))
        self.genomes_w_one = self.list_to_string(self.get_genomes(1))
        self.genomes_w_two = self.list_to_string(self.get_genomes(2))
        self.genomes_w_three = self.list_to_string(self.get_genomes(3))

        self.snp_total = self.get_snp_total()

        # A class variable for the assigned encoding to set later, externally from this class
        self.group = None

    #--------------------------------------------

    def list_to_string(self, l):

        if len(l) == 0:
            return "--"
        else:
            return ",".join(sorted(l))
    
    #--------------------------------------------

    def first_half(self):
        return self.snp_data[0:self.qindexes[0]]

    #--------------------------------------------

    def second_half(self):
        return self.snp_data[self.qindexes[0]:]
        
    #--------------------------------------------

    def ref_matches(self):
        """
        This method just returns the number of times the refbase
        matches all qbases.
        """
        return self.qbases.count( self.ref_base )

    #--------------------------------------------

    def set_informative(self):
        """
        Just going to use the global ones, twos and threes variables for the check.

        Simple check to find out if it's informative or not.
        """

        if self.ones == 1 or self.zeroes == 1:

            if self.twos == 0 or (self.twos == 1 and self.threes <= 1):

                return "NI"


        return "PI"

    #--------------------------------------------

    def get_genomes(self, num):
        """
        retreives all genomes in the self.query_genomes list that match up
        with the postion of the number in the pattern_list.
        """
        indexes = [i for i,val in enumerate(self.pattern_list[1:]) if val == num]

        return [self.query_genomes[i] for i in indexes]
        
    #--------------------------------------------

    def get_pattern(self):
        """
        performs the walk along the qbases and does the count if
        refbase bases match the qbase and resets if it doesn't.
        Only need to call this once on init since the SNP is immutable.
        """
        encoding = [0]
        seen = {}
        count = 0

        # set the first for refbase
        seen[self.ref_base] = 0
                
        for qb in self.qbases:

            if not seen.has_key(qb):
                count += 1
                seen[qb] = count


            encoding.append(seen[qb])


        return encoding


    #--------------------------------------------

    def get_num_refs(self, num):
        """
        Returns a list of the bases where all counts match 'num'
        """
        
        refs = []
        
        for indx, count in enumerate(self.pattern_list):

            if count == num:
                refs.append( self.qbases[indx] )

        return sorted(refs)

    #--------------------------------------------
    
    def get_snp_total(self):

        qbases = [q for q in self.qbases if not q in [self.ref_base, "No Hit", "indel"]]

        if len(qbases) == 0:

            return "--"

        else:

            qtotal = []
            
            # get rid of the stuff with a '/' in it.
            for q in qbases:

                if '/' in q:
                    
                    qtotal += [qb for qb in q.split('/') if not qb in [self.ref_base, "No Hit", "indel"]]

                else:

                    qtotal += q

            
            return "".join( list(set(qtotal)) )
        
    
#*** end SNP class *********************************************




#***** AlphaIncrementer class object ***********************************
class AlphaIncrementer:

    #--------------------------------------------
    def __init__(self, start=chr(ord('A') - 1)):

        self.current_char = 0


    #--------------------------------------------

    def get_next(self):
        """
        Increments the current character to the next one
        and returns the new current.
        """
        
        #nchar = chr(ord(self.current_char)+1)

        #self.current_char = nchar

        self.current_char += 1

        return str(self.current_char)
    
#****************************************************************




#****************************************************************
# 
#
#****************************************************************
class MoleculeDict:
    def __init__(self):
        self._dict = {}

    #--------------------------------------------
    def add(self, snp):

        if not self._dict.has_key(snp.pattern):

            self._dict[snp.pattern] = {}

        if not self._dict[snp.pattern].has_key(snp.molecule):

            self._dict[snp.pattern][snp.molecule] = []

            
        #self._dict[snp.pattern].append(snp)
        self._dict[snp.pattern][snp.molecule].append(snp)


    #--------------------------------------------
    def get(self, snp):

        reflist = []
        
        if self._dict.has_key(snp.pattern):

            pattern = self._dict[snp.pattern]

            
            
            for k in pattern.keys():

                refs = ",".join([str(s.ref_pos) for s in pattern[k]])

                molpos = "%s: %s" % (k,refs)

                reflist.append(molpos)

        

        return reflist


    #-----------------------------------------------
    def get_string(self, snp, unique=False):

        if len(self.get(snp)) > 0:

            items = []
            if unique:
                
                items = list(set(self.get(snp)))

            else:

                items = self.get(snp)

            return " | ".join(items)
            
        else:
            return "--"


#****************************************************************
# 
#
#****************************************************************
class GroupDict:
    def __init__(self):
        self._dict = {}


        self.molecules = {}

    #--------------------------------------------
    def add(self, snp):

        if not self._dict.has_key(snp.group):

            self._dict[snp.group] = {}
            self._dict[snp.group]['snp'] = snp

        if not self._dict[snp.group].has_key(snp.molecule):

            self._dict[snp.group][snp.molecule] = []
            
            

        if not self.molecules.has_key(snp.molecule):

            self.molecules[snp.molecule] = snp

            
        self._dict[snp.group][snp.molecule].append(snp)

    #--------------------------------------------

    def get_group_string(self, group, molecule_list):

        group_collection = self._dict[group]

        snp = group_collection['snp']



        molecule_refpos = [self.get_molecule_string(group, m) for m in molecule_list]
        
        line = "\t".join([str(group), str(snp.pattern), snp.info, snp.genomes_w_one, snp.genomes_w_two, snp.genomes_w_three, snp.genomes_w_zero] + molecule_refpos)
        
        return line 
    #--------------------------------------------        

    def get_molecule_string(self, group, molecule):

        if self._dict[group].has_key(molecule):

            snps = self._dict[group][molecule]

            return ",".join([str(s.ref_pos) for s in snps])

        else:

            return "--"
        
#****************************************************************

        
#-------------------------------------------------------------------------------

def load_table(table_file):
    """
    Loads the table file into a tuple that returns the header, qindexes
    and a list of all snp objects.
    """
    
    infile = open(table_file, 'r')

    table_data = infile.readlines()

    infile.close()

    # The header should always be the first line
    header = table_data[0].split('\t')

    # collect the qindexes from the header with this loop within a list
    qindexes = [indx for indx,colname in enumerate(header) if "qbase:" in colname]

    # collect the names of the genomes from the header
    query_genomes = [header[qi].replace("qbase:","").rstrip() for qi in qindexes]


    indexes = {}
    indexes['molecule'] = header.index('molecule')
    indexes['refpos'] = header.index('refpos')
    indexes['refbase'] = header.index('refbase')
    indexes['ref_codon'] = header.index('ref_codon')
    indexes['query_aa'] = header.index('query_aa')
    indexes['ref_aa'] = header.index('ref_aa')
    indexes['query_codon'] = header.index('query_codon')
    indexes['gene_name'] = header.index('gene_name')


    # put each line representing a SNP into a SNP object which was declared earlier
    # We exclude the first line since it's the header.
    snp_objects = []
    
    for snp_line in table_data[1:]: 

        snp_objects.append( SNP(qindexes, snp_line, query_genomes, indexes) )


    # returning everything in a tuple
    print "Loaded %d SNPs from file '%s' with %d qbases per SNP.\n" % (len(snp_objects), table_file, len(qindexes))
    
    return (header, qindexes, snp_objects)


#*******************************************************************************

def group_file(group_dict, outfile):


    with open(outfile, 'w') as of:
        
        molecules = sorted(group_dict.molecules.keys())

        header = ["Group", "Pattern", "Informative", "Genomes_w_1", "Genomes_w_2", "Genomes_w_3", ] + ["refs:%s" % m for m in molecules]

        
        groups = sorted(group_dict._dict.keys())


        of.write("\t".join(header)+"\n")


        for g in groups:

            line = group_dict.get_group_string(g, molecules)

            of.write(line+"\n")

    print "Write a group file to '%s'" % outfile

#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("snp_table", type=str,
                        help="The snp table to input")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="snp_genotype_out.txt")
    parser.add_argument("-t", "--translation_table", type=int, default=1,
                        help="The translation table to use for amino translations")
    parser.add_argument("-g", "--groupfile", type=str,
                        help="Output file indexed by pattern groups")


#    parser.add_argument("-m", "--min_table", action="store_true",
#                        help="Output a minimum table with only the patterns and codes.")

    args = parser.parse_args()

    output_file = args.outfile

    # open the snp table and load it into some data types.
    
    header, qindexes, snp_objects = load_table(args.snp_table)

    # Here I'm going to encode the snp objects in a hash
    # using the alpha incrementer.
    alphainc = AlphaIncrementer()
    encoding_dict = {}
    molecule_dict = MoleculeDict()
    group_dict = GroupDict()

    # loop through snp objects and set each one
    for snp in snp_objects:

        # if this is not given a key, then we add it and set it.
        if not encoding_dict.has_key(snp.pattern):

            encoding_dict[snp.pattern] = alphainc.get_next()

        snp.group = encoding_dict[snp.pattern]


        # add this snp info to the molecule and refpos dicts
        molecule_dict.add(snp)
        group_dict.add(snp)
        
    #------- Finished loop for encoding ---------------------------------------


    #
    # From here we just output with the new columns, same loop order as before.
    #
    with open(output_file, 'w') as of:

        newcols_start = qindexes[0] 
        newcols_end = qindexes[0] 

        # First write the header.
        of.write( "\t".join(header[0:newcols_start] + ["snp_total", "Pattern", "Group", "Informative", "Genomes_w_1", "Genomes_w_2", "Genomes_w_3", "Genomes_w_zero", "Reference Positions"] + header[newcols_end:]) )

        for snp in snp_objects:

            line = "\t".join( snp.first_half() + [ snp.snp_total, snp.pattern, snp.group, snp.info, snp.genomes_w_one, snp.genomes_w_two, snp.genomes_w_three, snp.genomes_w_zero, molecule_dict.get_string(snp)] + snp.second_half() ) 

            of.write(line)

    print "Output %d SNPs to file %s" % (len(snp_objects), output_file)


    # If we want to write a group file out we do it here.
    if not args.groupfile is None:
        
        group_file(group_dict,args.groupfile)

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
