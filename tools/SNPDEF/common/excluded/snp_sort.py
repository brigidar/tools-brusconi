#!/usr/bin/env python

import argparse, os, sys
import pdb


import operator

from Bio import SeqIO



#-------------------------------------------------------------------------------
#
# Takes a snp in the form of a string and places it into a class
# data structure.
#-------------------------------------------------------------------------------
class SNP:
    
    def __init__(self, molecule, pos, line):
        
        self.molecule = molecule
        self.pos = int(pos)
        self.line = line

#-------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------
class SNPTable:
    
    def __init__(self, table_file, out_file, exclude_file=None, exclude_list=None):
        
        
        self.header = ""
        
        # This will be a dict with a list of tuples in the form exclude_list[molecule] = [(start, stop)...]
        self.exclude_list = {}
        
        self.snps = []
        
        self.read_table_file(table_file)
        
        if not exclude_file is None:
            self.read_exclude_file(exclude_file)
        
        
        print "Filtering %d snps" % (len(self.snps))
        
        self.filtered_snps = []
        
        for snp in self.snps:
            
            if self.in_exclude(snp):
                continue
            else:
                self.filtered_snps.append(snp)
        
        
        print "Sorting %d snps" % (len(self.filtered_snps))
        
        self.filtered_snps.sort(key=lambda x: x.pos)
        
        
        self.write_file(out_file)

#-------------------------------------------------------------------------------

    def read_table_file(self, tf):
    
        with open(tf, 'rU') as input_handle:
        
            for line in input_handle:
            
                parts = line.rstrip().split('\t')
                
                if parts[0] == "molecule" and parts[1] == "refpos":
                    
                    self.header = line
            
                else:
                    
                    self.snps.append( SNP(parts[0], int(parts[1]), line) )

#-------------------------------------------------------------------------------

    def in_exclude(self, snp):
    
        try:
        
            for s in self.exclude_list[snp.molecule]:
            
                if s[0] < snp.pos < s[1]:
                    return True
        except KeyError:
            return False
    
        return False
    
    #-------------------------------------------------------------------------------
    
    def read_exclude_file(self, exclude_file):
        
        with open(exclude_file, 'rU') as input_handle:
            
                for l in input_handle.readlines():
                
                    parts = l.rstrip().split('\t')
                
                    if parts[0] == "molecule":
                        continue
            
                    else:
                    
                        self._add_to_exclude_list(parts[0], int(parts[1]), int(parts[2]))

#-------------------------------------------------------------------------------

    def _add_to_exclude_list(self, molecule, start, stop):
    
        if self.exclude_list.has_key(molecule):
        
            self.exclude_list[molecule].append( (start, stop) )
        
        else:
            
            self.exclude_list[molecule] = [ (start, stop) ]

#-------------------------------------------------------------------------------

    def write_file(self, out_file):
    
        with open(out_file, 'w') as output_handle:
        
            output_handle.write(self.header)
            
            for fs in self.filtered_snps:
                
                output_handle.write(fs.line)
        
                print "Done writing the output to %s" % out_file

#-------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------

def __main__():
    
    
    #Parse Command Line
    parser = argparse.ArgumentParser()
    
    
    # Might have to change these
    parser.add_argument('-s', '--snp_table', required=True, help="a snp table file", metavar="snp_panel_file")
    parser.add_argument('-f', '--exclude_file', required=True, help="a template file", metavar="exclude_file")
    parser.add_argument('-o', '--out', help="Output file for the merged table", metavar="sorted_snp_output_file")
    
    
    
    
    args = parser.parse_args()
    
    #pdb.set_trace()
    snp_table = SNPTable(table_file=args.snp_table,out_file=args.out,exclude_file=args.exclude_file)
                         
    print "Done filtering %s to %s" % (args.snp_table, args.out)


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
