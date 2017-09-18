#!/usr/bin/env python

__version__="1.15.0"

"""
Name: snp_verify.py

Description : Script to parse BLAST results and annotate verified SNPs.

Author: Mando Rodriguez with consulting from Brigida Rusconi

Based on the snp_verify.pl script in Ergartis written by Sonia Agrawl
for IGS. Program features options for multithreading and multiprocesing.


"""


import argparse, os, sys
import pdb
import operator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import ExtendedIUPACDNA
from Bio.SeqFeature import FeatureLocation

import multiprocessing
import logging
from multiprocessing.dummy import Pool as ThreadPool


# these import are for timekeeping
import timeit


"""

seqrecord file:
    http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html

the hsp class:
    http://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp.HSP-class.html

possible multiprocessing:
    http://pymotw.com/2/multiprocessing/basics.html


Some code for different protein translations

http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2

"""

#-------------------------------------------------------------------------------
# function just lets you translate a codon based on a specific amino table
#
# May need to make the type selectable for the Seq, details here:
#     http://biopython.org/wiki/Seq
#-------------------------------------------------------------------------------
def translate_codon(codon, table=1):

    if codon == "--":
        return "--"
    elif codon == "indel":
        return "indel"
    else:
        return str(Seq(str(codon), ExtendedIUPACDNA()).translate(table=table))

#-------------------------------------------------------------------------------
# Function gives whether a nucleotide change is a transition or a transversion
#
#-------------------------------------------------------------------------------
def get_trans_state(base, hit):

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

        return "Error [base:%s, hit:%s]" % (base, hit)

#-------------------------------------------------------------------------------
# This function computes a set of contigs and sets it to the snp
#-------------------------------------------------------------------------------

def compute_snp(sd):
    """
    This is the function that is passed to the threaded and multiproc
    branches to produce snps.
    """
    try:
        return SNP(sd[0], sd[1], sd[2], sd[3], sd[4], sd[5], sd[6], sd[7])
    except Exception, e:
        print "Unable to create a snp entry for molecule '%s' at snp position %d: %s" % (sd[0], sd[1], e)
        #logger.exception("Unable to create a snp entry for molecule '%s' at snp position %d: %s" % (sd[0], sd[1], e))
        raise Exception(e)


#-------------------------------------------------------------------------------


class Hit:
    """
    This class collects the base nucleotide, codon, offset, and blast hit into a class
    for storing and outputting. All data from a hit collected by a SNP object is contained in this object.
    """
    def __init__(self, ref_base, ref_amino, base, codon, offset, blast_hit, translation_table=1):
        self.base = base
        self.is_snp = str(ref_base) != base
        self.offset = offset
        self.blast_hit = blast_hit
        self.codon = codon
        self.aa = None
        
        if not codon is None:
            self.aa = translate_codon(codon,translation_table)
            
        self.trans = get_trans_state(str(ref_base), str(base)) 

        self.stat = "NA"
        if ref_amino == self.aa:
            self.stat = "SYN"
        elif self.aa is None:
            self.stat = None
        elif self.aa == "--":
            self.stat = "--"
        elif self.aa == "indel":
            self.stat = "indel"
        else:
            self.stat = "NSYN"
    
    def length(self):
        return str(self.blast_hit.hit_span)
    
    def __str__(self):
        return self.base
    
#-------------------------------------------------------------------------------
# 
# 
#-------------------------------------------------------------------------------
class SNP:

    """
    SNP class will take all of the relevant data to search through and produce a verified snp upon initialization.
    
    """
    
    def __init__(self, molecule, pos, contig_files, blast_data, genbank_refs,
                 flanking_bases=None, blast_threshold=0.0, blast_length_threshold=0,
                 amino_table=1):

        # Here I parse apart the molecule section given from the snp panel in the event that
        # there are pipes or dots in the identifier.

        self.molecule = molecule
        self.molecule_name = molecule
        self.molecule_parts = [mp for mp in molecule.split('|') if mp != '']


        # Here I look for anything with a dot in them, and from there I will append an identifier
        # without the dot into molecule parts.
        for mp in self.molecule_parts:

            if '.' in mp:
                try:
                    self.molecule_parts.append( mp.split('.')[0] )
                except Exception,e:
                    print "Error taking prefix from dot in molecule identifier %s: %s" % (mp,e)

        # Here's the whole thing to check for a match using only the 'in' operator.
        self.molecule_parts.append(molecule)
        
        self.pos = pos
        self.name = None
        self.id = None
        self.tag = "%s_SNP_%d" % (molecule, pos)

        # Gene data for this SNP
        self.ref_codon = None
        self.ref_amino_acid = None
        self.ref_base = None
        self.ref_pos = pos - 1 
        
        self.gene = None
        self.gene_name = "intergenic"
        self.gene_start = None
        self.gene_end = None
        self.gene_length = None
        self.snps_per_gene = None
        self.pos_in_gene = None
        self.pos_in_codon = None

        self.genome_end = None

        self.codon_start = None

        self.product = "intergenic"
        self.note = "No Note"
        self.strand = None

        self.base_hits = []
        self.total_hits = 0
        
        # these are for sections that coalese column data
        self.all_codons = []
        self.all_aas = []
        self.all_trans = []
        self.all_stats = []


        # Here we add all of the gene data to the SNP
        self.collect_gene_info(pos, genbank_refs)

        # here we translate the ref codon that the previous method call collected
        self.ref_amino_acid = translate_codon(self.ref_codon, amino_table)

        # Now we go through all of the queries to collect the data in the order of the contig_files so they are
        # in the table order.
        self.query_contigs = [QueryContig(self, cf, blast_data, flanking_bases, blast_threshold, blast_length_threshold, self.ref_base) for cf in contig_files]


        # Here we collect all of the hits we gathered in the contig files to coalese them into one column
        # for certain data.
        for qc in self.query_contigs:
            self.base_hits += qc.base_hits

        # in place sort
        self.base_hits.sort(key=lambda x : x.base)


        # here we collect the grouped data for the columns that need it
        for bh in self.base_hits:

            if bh.is_snp:
                self.total_hits += 1

                if not bh.codon is None:
                    self.all_codons.append(str(bh.codon))
                    
                self.all_trans.append(bh.trans)
                self.all_stats.append(bh.stat)
                self.all_aas.append(bh.aa)


    #-------------------------------------------------------------------------------

    def molecule_name(self):

        # If we have pipes in the id then we keep the third entry for molecule
        if len(self.molecule_parts) >= 3:
            try:
                return self.molecule_parts[2]
            except Exception:
                return self.molecule

    #-------------------------------------------------------------------------------

    def get_gene_start(self):
        if self.gene_name == "intergenic":
            return "NA"
        else:
            return str(self.gene_start+1)

    #-------------------------------------------------------------------------------

    def get_gene_end(self):
        if self.gene_name == "intergenic":
            return "NA"
        else:
            return str(self.gene_end+1)

    #-------------------------------------------------------------------------------

    def get_pos_in_gene(self):
        if self.gene_name == "intergenic":
            return "--"
        else:
            return str(self.pos_in_gene+1)

    #--------------------------------------------------------------------------------

    def get_stats(self):
        """
        Get's all stat strings from the indivdual query entries and outputs them into a string.
        """

        if self.gene_name == 'intergenic':

            return 'intergenic'

        if len(self.all_stats) == 0:

            return "--"

        else:
            
            return '/'.join(sorted(self.all_stats))



        
    #-------------------------------------------------------------------------------

    def get_codons(self):
        """
        Returns all collected codons in the SNPs
        """

        if len(self.all_codons) == 0:

            return "--"

        else:
            
            return '/'.join(sorted(self.all_codons))

    #-------------------------------------------------------------------------------

    def get_aas(self):
        """
        Returns a list of all translated amino acids collected in the SNPs
        """

        if self.gene_name == "intergenic":

            return "intergenic"
        
        if len(self.all_aas) == 0:

            return "--"

        else:
            
            return '/'.join(sorted(self.all_aas))

    #-------------------------------------------------------------------------------

    def get_trans(self):
        """
        Returns all the trans states collected in the verified SNPs
        """

        if len(self.all_trans) == 0:

            return "--"

        else:
            
            return '/'.join(sorted(self.all_trans))


    #-------------------------------------------------------------------------------

    def collect_gene_info(self, snp_pos, genbank_refs):
        """
        This method will take the molecule and snp position and
        get all of the gene info and add it to the snp object.
        """
        
        source = self.get_source_from_parent(self, genbank_refs)

        self.name = source.name
        self.id = source.id
        self.strand = 1 # assumed upstream unless it finds a gene

        # Save the genome end of the source.
        self.genome_end = len(source)-1
        
        pos = self.ref_pos

        gene_hits = self.find_gene_hits(pos, genbank_refs)

        if len(gene_hits) > 0:

            # just going to use the first gene. Shouldn't be any overlaps with contigs.
            gene_tuple = gene_hits[0]

            self.add_gene_info(self, gene_tuple)


        #-------- finished getting gene info ------------------------------------
        
        try:

            self._get_codon_and_base(source.seq)

        except Exception, e:
            raise Exception("Error getting reference (base,codon,aa) for %s at position %s: %s" % (self.molecule, snp_pos, e))


    #-------------------------------------------------------------------------------

    def add_gene_info(self, snp, gene_tuple):
        """
        If the gene is a hit, then we add some identifiers to the snp object
    
        This is handled here instead of inside the snp object because we need to fetch
        data from the parsed sequences
        """
        
        snp.gene = gene = gene_tuple[0]
        snp.name = gene_tuple[1]
        snp.id = gene_tuple[2]
        
        try:
            
            locus_tags = gene.qualifiers['locus_tag']

            gene_name = gene.qualifiers['locus_tag'][0]

        except Exception, e:
            raise Exception("No genes (locus_tag) for gene hit on molecule '%s' at (%d,%d)" % (snp.molecule_name(), int(gene.location.start)+1, gene.location.end))
        
        snp.strand = gene.strand
        snp.gene_length = str(len(gene.location))
        snp.gene_name = gene_name

        pos = snp.ref_pos

        # accomodate for strands
        #
        #  The location end must have 1 subtracted to have it
        #  point to the exact position.
        #
        if snp.strand == 1:

            snp.gene_start = int(gene.location.start)
            snp.gene_end = int(gene.location.end)
            snp.pos_in_gene = (pos - snp.gene_start)
            

        elif snp.strand == -1:

            snp.gene_start = int(gene.location.end)-1
            snp.gene_end = int(gene.location.start)
            snp.pos_in_gene = (snp.gene_start - pos)


        snp.pos_in_codon = (snp.pos_in_gene+1) % 3

        # zero means it's at the third base of the codon since it divides evenly
        # into three.
        if snp.pos_in_codon == 0:

            snp.pos_in_codon = 2

        else:

            snp.pos_in_codon = snp.pos_in_codon - 1

        try:
            snp.product = ','.join(gene.qualifiers['product'])
        except KeyError:
            snp.product = "intergenic"

        try:
            snp.note = ','.join(gene.qualifiers['note'])
        except KeyError:
            snp.note = "No Note"


    #-----------------------------------------------------------------------------
 
    def _get_codon_and_base(self, seq):
        """
        gets the codon at that position and tacks on 'N' for any blank entries. It updates class
        members interally.

        It doesn't calc the amino since that is done after the other data is collect by way of the codon.
        """
        
        pos = self.ref_pos

        if self.gene_name == "intergenic":

            #codon =  FeatureLocation(pos-snp.pos_in_codon, pos-snp.pos_in_codon+3, 1).extract(seq)

            self.ref_base = str(FeatureLocation(pos, pos+1, strand=self.strand).extract(seq))

            self.ref_codon = "--"

            self.codon_start = -1
        
        else:

            # The gene location object already has the positions and the strand variables, so if it needs
            # to be complemented so the gene_seq should be complemented. 
            gene_seq = self.gene.location.extract(seq)

            start = 0

            if self.ref_pos == self.gene_start:

                start = 0

            else:
                
                start = self.pos_in_gene - self.pos_in_codon

                if start < 0:

                    raise Exception("Error getting start position for %s: pos_in_gene(%s), pos_in_codon(%s), calculated start is %s" % (self.tag, self.pos_in_gene, self.pos_in_codon, start))

            
            self.ref_codon = FeatureLocation(start, start+3).extract( gene_seq )

            self.ref_base = str(FeatureLocation(self.pos_in_gene, self.pos_in_gene+1).extract( gene_seq ))

            self.codon_start = start

    #-------------------------------------------------------------------------------
    
    def get_source_from_parent(self, snp, genbank_refs):
        """
        Retrieves the source sequence from the genbanks
        """
        
        # Here we go through the genbanks to get the intergenic sequence
        for gblist in genbank_refs:

            for g in gblist:

                if g.name in snp.molecule_parts or g.id in snp.molecule_parts:

                    # returns the first matching genbank record
                    return g

        raise Exception("Can't extract source sequence at pos %d from molecule %s, no matching genbanks" % (snp.pos, snp.molecule))

    #-------------------------------------------------------------------------------

    def get_gbank_hits(self,genbank_recs,molecule_id):

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
    # Returns a list of genes and gene locations that the refpos lies within
    #
    #-------------------------------------------------------------------------------
    def find_gene_hits(self, refpos, genbank_refs):
        """
        retrieves all of the gene hits on a molecule and reference position from the genbank
        features. If it's on the negative strand the end must be inclusive while the start is not.
        """
        gene_hits = []
            
        for gblist in genbank_refs:

            for gb in gblist:

                #if molecule in [gb.name,gb.id]:
                if gb.name in self.molecule_parts or gb.id in self.molecule_parts:

                    # need to save the one that matches into this field for output
                    self.molecule_name = gb.name if gb.name in self.molecule_parts else gb.id
                    
                    for feature in gb.features:

                        location = feature.location

                        if feature.type == "CDS" :

                            # on the positive strand
                            if location.strand == 1 and (location.start <= refpos < location.end):

                                    # here we append the whole gene with the hit
                                    gene_hits.append((feature, gb.name, gb.id))

                            # on the negative strand
                            elif location.strand == -1 and (location.start < refpos <= location.end-1):

                                    gene_hits.append((feature, gb.name, gb.id))                                

                               
        return gene_hits

#-------------------------------------------------------------------------------





#**************************************************************************************

class QueryEntry:
    """
    This class is a 'child' class to class QueryContig.
    
    For each entry it keeps track of the file the query is in
    as well as the name/id of the genome that it is collecting data for.
    The hits are stored in the dict self.base_hits
    
    """
    def __init__(self, snp, fasta_record, blast_hits, flanking_bases, translation_table=1):

        self.snp = snp
        self.fasta = fasta_record # just to keep the data around
        self.name = fasta_record.name
        self.id = fasta_record.id
        self.blast_hits = blast_hits # these are only matching now.

        self.flanking_bases = flanking_bases
        self.translation_table = translation_table
        
        self.base_hits = []
        self.shortest_blast = "NA"
        self.longest_blast = "NA"

        self.collect_base_hits()

    #-------------------------------------------------------------------------------

    def calc_offset(self, pos, hit):
        """
        Tries to produce the offset needed to extract the nucleotide
        at the position indicated in 'pos'

        Will use data from the blast hit to pull this out since we have
        all the proper offsets relative to the contig saved in the query
        section of the hit.

        """
        
        nuc_offset = pos if pos < self.flanking_bases else self.flanking_bases

        # if the query start on the hit is greater than zero
        # then that means the query.seq is shorter in the front
        # by that amount and we need to reduce the offset to get the
        # correct position.
        
        nuc_offset -= hit.query_start
        
        if str(hit.query.seq).count('-') > 0:

            indel_indexes = [i for i,k in enumerate(str(hit.query.seq)) if k == '-']

            for i in indel_indexes:

                if i < nuc_offset:
                    nuc_offset += 1
                            
            
        return nuc_offset
    
    #-------------------------------------------------------------------------------

    def calc_offset_from_left(self, pos, hit):
        """
        This is the older less streamlined version for getting the base offset
        but since it was always correct I'm keeping it around to compare in case
        any issues come up with the newer offset calculation.

        """

        query_start = pos-self.flanking_bases 
        query_end = query_start + (self.flanking_bases*2 + 1)

        hit_start = query_start + hit.query_start
        hit_end = hit_start + hit.query_span

        offset = pos-hit_start

        if str(hit.query.seq).count('-') > 0:

            indel_indexes = [i for i,k in enumerate(str(hit.query.seq)) if k == '-']

            for i in indel_indexes:

                if i < offset:
                    offset += 1
                            
            
        return offset

    #-------------------------------------------------------------------------------

    def collect_base_hits(self):
        """
        this method will use the flanking base and match it up against the snp and the blast
        hits.
        """
        pos = self.snp.ref_pos

        hit_base = None
        query_base = None
        
        for bh in self.blast_hits:

            offset = -1

            offset = self.calc_offset(pos, bh)
            
            if pos < self.flanking_bases:
                
                # this is just a message to notify the user about a nucleotide close to the beginning.
                print "SNP %s is near the beginning of the genome, at offset %d" % (self.snp.tag, offset)
                

            try:

                if offset < 0:
                    raise IndexError("Calculated No Hit on query positions (%s,%s)" % (bh.query_start, bh.query_end))


                hit_base = None
                query_base = None
                    
                # this is just here to trigger an exception if out of range
                bh.hit.seq[offset]

                if self.snp.strand == 1:

                    hit_base = bh.hit.seq[offset:offset+1]
                    query_base = bh.query.seq[offset:offset+1]
                    
                elif self.snp.strand == -1:

                    hit_base = bh.hit.seq[offset:offset+1].complement()
                    query_base = bh.query.seq[offset:offset+1].complement()
                    

                # here we collect the shortest blast hit
                if self.shortest_blast == "NA" or bh.hit_span < self.shortest_blast:

                    self.shortest_blast = bh.hit_span

                # this gets the highest blast
                if self.longest_blast == "NA" or bh.hit_span > self.longest_blast:

                    self.longest_blast = bh.hit_span


                # Here we set the hit base
                if str(hit_base) == "-":
                    hit_base = "indel"
                else:
                    hit_base = str(hit_base).upper()


                try:
                    
                    qcodon = self.get_qcodon(hit_base)

                    self.base_hits.append( Hit(self.snp.ref_base, self.snp.ref_amino_acid, hit_base,
                                               qcodon, offset, bh, self.translation_table) )

                except Exception,e:
                    raise Exception("Error obtaining qcodon/qaa for %s and FASTA record '%s': %s" % (self.snp.tag,self.id,e))
                
            except Exception, e:

                print "Can't get base for query '%s' from blast '%s:%s' with seq length '%d'and offset '%s': %s" % (self.name, bh.query_id, bh.hit_id, len(bh.hit.seq),offset, e)

        
    #-------------------------------------------------------------------------------

    def get_qcodon(self, qbase):
        """
        Fetches the query codon to the query contig. If a
        base is the same as the reference it will not be added to the hits.
        
        """
        
        if self.snp.gene_name == "intergenic":

            # if it's intergenic we do nothing.
            return
        
        if str(qbase) == "indel":
            

            return "indel"
            
        else:
            pic = self.snp.pos_in_codon
            
            ref_codon = self.snp.ref_codon # already complemented if necessary
            
            ql = list(str(ref_codon))

            try:
                ql[pic] = str(qbase) #should never fail
            except Exception:
                
                raise Exception("Can't get query codon for ref_codon %s and query base %s with codon pos %d" % (ref_codon, qbase, pic))

            # collect the qcodons only if it differs from the reference
            
            qcodon = Seq(''.join(ql), ref_codon.alphabet)

            # Don't really need to do this here since we do it at the SNPVerify
            # level, but i'm leaving it here for the error to catch. We don't need it
            # since we just calculate all aminos from the total list of codons at the very
            # end of the run for the snps.
            qaa = qcodon.translate()

            return qcodon


    
#************************ End QueryEntry ***************************************




#********************************************************************************
# This is the data loaded from a file full of contigs. It is meant to be immutable
# and processes will simply search through it. 
#
#********************************************************************************
class ContigFile:
    """
    This is the data loaded from a file full of contigs. The object is meant to be immutable
    and processes will simply search through it. 
    """
    def __init__(self, filename):
        
        self.base_filename = ""
        self.meta_filename = ""
        self.contig_records = []

        lineparts = filename.split(':')

        if len(lineparts) == 2:

            self.base_filename = lineparts[0]
            self.meta_filename = lineparts[1].replace(' ','_')

        elif len(lineparts) == 1:

            self.base_filename = lineparts[0]

            try:
                self.meta_filename = os.path.splitext(os.path.basename(lineparts[0]))[0]
            except Exception:
                self.meta_filename = lineparts[0]

        else:
                    
            raise Exception("Can't parse '%s', not a valid query contig FASTA file" % filename)


        # here we load up the fasta records

        fasta_seqgen = SeqIO.parse(open(self.base_filename), 'fasta')

        for f in fasta_seqgen:

            self.contig_records.append(f)

        print "Parsed query FASTA file: %s with %d fasta records" % (filename, len(self.contig_records))


#*********************************************************************************

class PickleableQueryResult:
    """
    This is a generic query result, only it's pickleable since it
    removes the function arguments and lambda functions. All we need is the
    'hits' section.
    """
    def __init__(self, qrid, hits):
        self.id = qrid
        self.hits = hits


#**************************** class QueryContig **********************************
class QueryContig:
    """
    This class is a 'child' class to class SNP.

    It is composed of a listof QueryEntry classes that compile the final data
    for self.base_hits and self.base_count respectively.

    Each QueryContig takes a QueryFile object as a constructor argument to keep
    track of the file they are in.
    
    This does the legwork for comparing the queries to the snp. It's a direct mapping; 
    one query contig (fasta) to one ContigFile. 
    """
    def __init__(self, snp=None, contig_file=None, blast_data={}, flanking_bases=None,
                 blast_threshold=0.0, blast_length_threshold=0, ref_base=None):

        #-- now plug in our values --
        self.snp = snp
        self.query_id = "%s_%s_SUBSEQ" % (self.snp.molecule, self.snp.pos) 
        self.base_filename = contig_file.base_filename
        self.meta_filename = contig_file.meta_filename
        self.contig_records = contig_file.contig_records
        self.blast_threshold = blast_threshold
        self.blast_length_threshold = blast_length_threshold
        self.flanking_bases = flanking_bases
        self.ref_base = ref_base

        self.longest_blast = 0

        self.matching_queries = self.get_matching_query_results(blast_data)

        # a list of collected querie entries, will not only be relagated
        # to hits.
        self.queries = []

        # collect the query entries
        self.get_query_entries()

        # just to cache the base hits so we can return this instead.
        self.base_hits = []

        # Here we collect all hits from the query entries. The base hits are the only thing we need
        # from the entries. 
        for qe in self.queries:
            
            self.base_hits += qe.base_hits

            if qe.longest_blast > self.longest_blast:
                self.longest_blast = qe.longest_blast

        # in place sort
        self.base_hits.sort(key=lambda x : x.base)

        self.valid_hits = [bh for bh in self.base_hits if bh.is_snp]

        self.total_hits = len(self.valid_hits)

    #-------------------------------------------------------------------------------

    def get_base_nucs(self):
        """
        Returns a string conaining the nucleotides from the hits separated by a slash
        """

        base_hits = [bh.base for bh in self.base_hits]

        if len(base_hits) == 0:

            return "No Hit"

        else:

            return '/'.join(base_hits)

    #-------------------------------------------------------------------------------

    def get_blast_lengths(self):
        """
        
        """

        blast_lengths = [bh.length() for bh in self.base_hits]

        if len(blast_lengths) == 0:

            return "No Hit"

        else:

            return '/'.join(blast_lengths)
        

    #-------------------------------------------------------------------------------

    def get_query_entries(self):
        """
        This method loops through all contigs and gets blast hits. If there are
        hits then it will place the snp, contig fasta record, and blast hits
        into a query entry.
        """
        
        for cr in self.contig_records:

            bhits = self.get_blast_hits(cr)

            if len(bhits) == 0:
                # if no hits then we just keep going and don't save anything.
                continue

            self.queries.append( QueryEntry(self.snp, cr, bhits, self.flanking_bases) )

    #-------------------------------------------------------------------------------

    def get_matching_query_results(self, blast_data):
        """
        compiles a list of query results from the blast records that match
        our SNP. Simple hash to the value.
        """
        query_results = []
        
        try:

            query_results = blast_data[self.query_id]

        except KeyError:
            pass

        return query_results
    
    
    #-------------------------------------------------------------------------------


    def get_blast_hits(self, contig):
        """
        Goes through all of the query results and looks for matches for the hits
        for each contig in this set. On a match it creates a QueryEntry and adds it to
        the list. Hits can be filtered by the bitscore and by the hit_span length
        """
        blast_hits = []

        for query_result in self.matching_queries:

            # look through a list of query results that match
            # the query name to the query result id
            for hit in [h for h in query_result.hits if h.id == contig.name]:

                for hsp in hit.hsps:

                    if hsp.query_id == self.query_id and hsp.hit_id == contig.id:

                        if hsp.bitscore >= self.blast_threshold and hsp.hit_span >= self.blast_length_threshold:

                            blast_hits.append(hsp)


        return blast_hits
    
    #-------------------------------------------------------------------------------

    def get_hit_ids(self):
        """
        Returns a list of the blast hit ids collected in the blasts.
        """
        hit_ids = []
        
        for qe in self.queries:

            for bh in qe.binfo:
                
                hit_ids.append(bh.hit_id)

        return hit_ids
    
#**************************** End class QueryContig ********************************
    

#**************************** The overall snp verify class *************************        
class SNPVerify:

    """
    may end up using a pool of workers for the contig processing

    https://docs.python.org/2/library/multiprocessing.html#using-a-pool-of-workers
    """
    #-------------------------------------------------------------------------------

    def __init__(self, flanking_bases=20, blast_threshold=None, blast_length_threshold=0,
                 amino_table=1,
                 multiproc=None, threads=None, 
                 blast_list_file=None, genbank_list_file=None, snp_list_file=None, query_contig_list_file=None, 
                 blast_list=[], genbank_list=[], snp_list=[], query_contig_list=[]):


        self.flanking_bases = int(flanking_bases)

        print "Using flanking base %d" % self.flanking_bases

        self.amino_table = int(amino_table)

        print "Using amino table %d" % self.amino_table
        
        self.merged_table = []
        self.genbank_refs = []

        self.snps_per_gene = {}
        self.syn_nsyn_per_gene = {}

        # This holds the entire table before printing it out
        self.verified_snps = []


        # store filenames
        # Possibly check to see if the file exists and is valid
        self.snp_list_file = snp_list_file
        self.query_contig_list_file = query_contig_list_file
        self.blast_list_file = blast_list_file
        self.genbank_list_file = genbank_list_file

        # here are lists read in from the command line
        self.blast_list = blast_list
        self.genbank_list = genbank_list
        self.snp_list = snp_list
        self.query_contig_list = query_contig_list

        try:
            self.blast_threshold = float(blast_threshold)
        except Exception:
            self.blast_threshold = 0.0

        if self.blast_threshold > 0.0:
            print "Using blast bitscore threshold of %d" % self.blast_threshold


        self.blast_length_threshold = blast_length_threshold

        if blast_length_threshold < 0 or blast_length_threshold > (self.flanking_bases*2 + 1):
            self.blast_length_threshold = (self.flanking_bases*2 + 1)
        
        print "Using blast length threshold %d" % self.blast_length_threshold


        #-- just some book keeping to see how many hits are filtered --
        self.num_blast_hits_loaded = 0
        self.num_filtered_blast_hits = 0
            
        self.multiproc = multiproc

        self.threads = None
        if not threads is None and threads != 0:
            self.threads = threads
            
        #if not threads is None and threads < 1:
            
            #print "Invalid value for threads (%d), must be a positive number" % threads
            #sys.exit()

        #else:
        #    self.threads = threads

        #------------------- high memory ----------------------------------------------
        # Here's some data that will take up a lot of memory when it loads
        #------------------------------------------------------------------------------
        self.snp_positions = {}
        self.blast_data = {}

        # a list for the query contigs from files
        self.contigs = []
        #------------------------------------------------------------------------------
        
        # build the table
        self.build_table()


    #-------------------------------------------------------------------------------

    def table_header(self):

        part1 = ["molecule", "refpos"] 

        part2 = ["syn/nsyn/intergenic"] # Now this will just be one column

        part3 = ["refbase"]

        part4 = ["qbase:%s" % c.meta_filename for c in self.contigs]

        part5 = ["gene_name",
         "gene_start", "gene_end", "gene_length", "snps_per_gene",
         "pos_in_gene","ref_codon", "ref_aa"
         ]
        
        part6 = ["query_codon"]

        part7 = ["query_aa"]

        part8 = ["transition/transversion"]

        part9 = ["snps/gene_length"]

        part10 = ["dn/ds"]

        part11 = ["num_hits:%s" % c.meta_filename for c in self.contigs]

        part12 = ["maxlen:%s" % c.meta_filename for c in self.contigs]

        part13 = ["blengths:%s" % c.meta_filename for c in self.contigs]

        part14 = ["product"]

        snp_list = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9 + part10 + part11 + part12 + part13 + part14
        
        return "\t".join(snp_list)


    #-------------------------------------------------------------------------------
        
    def snp_to_string(self, snp):

        part1 = [str(snp.molecule_name), str(snp.pos)] 

        part2 = [ snp.get_stats() ]

        part3 = [str(snp.ref_base)]

        part4 = [qc.get_base_nucs() for qc in snp.query_contigs]
            
        part5 = [str(snp.gene_name), str(snp.get_gene_start()), str(snp.get_gene_end()),
                 str(snp.gene_length), self.get_snps_per_gene(snp.molecule_name, snp.gene_name),
                 str(snp.get_pos_in_gene()),str(snp.ref_codon), str(snp.ref_amino_acid)]
        
        part6 = [snp.get_codons() ]

        part7 = [ snp.get_aas() ]

        part8 = [ snp.get_trans() ]

        part9 = [self.get_snp_div_gene_length(snp)]
        part10 = [self.get_syn_to_nsyn_ratio(snp)]

        part11 = [str(qc.total_hits) for qc in snp.query_contigs]

        part12 = [str(qc.longest_blast) for qc in snp.query_contigs]

        part13 = [qc.get_blast_lengths() for qc in snp.query_contigs]
        
        part14 = [str(snp.product)]

        snp_list = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9 + part10 + part11 + part12 + part13 + part14

        return "\t".join(snp_list)


    #-------------------------------------------------------------------------------
    # 
    #
    #
    #-------------------------------------------------------------------------------
    def build_table(self):

        # Here we parse each of the files, they get loaded into class
        # variables.
        #
        # First we load up the genbank files, snp positions, blast records, and query genomes.
        self.parse_genbank_list()
        self.parse_snp_positions()
        self.parse_blast_data()
        self.parse_query_genomes()

        
        print "Attempting to buld merged table\n"

        if not self.multiproc is None:

            # This is just to set the logging option
            #logger=multiprocessing.log_to_stderr(logging.DEBUG)
            #logger.setLevel(20)
            
            #--------------------------------------------------------
            # this uses a process pool to compute the jobs
            #--------------------------------------------------------


            # here we determine how many cpus to use
            if self.multiproc < 0 or self.multiproc > multiprocessing.cpu_count():
                self.multiproc = multiprocessing.cpu_count()
            
            print "Using %d of %d available cpus" % (self.multiproc, multiprocessing.cpu_count())


            inputs = []
            for molecule, snp_positions in self.snp_positions.iteritems():

                print "Setting up %d SNP jobs for molecule '%s'" % (len(snp_positions),molecule)

                for snp_pos in snp_positions:
                
                        inputs.append((molecule, snp_pos,
                                       self.contigs, self.blast_data, self.genbank_refs,
                                       self.flanking_bases, self.blast_threshold, self.blast_length_threshold))


            print "Running %d individual input jobs for processing on %d processes..." % (len(inputs), self.multiproc)
            

            print "Running multiprocessing mode..."
            pool = multiprocessing.Pool(self.multiproc)
            pool_output = pool.map(compute_snp, inputs)
            pool.close()
            pool.join()
            
            print "Collecting results..."

            # clearing out this memory
            inputs = None
            
            self.verified_snps = [r for r in pool_output]

                    
        elif not self.threads is None:

            #--------------------------------------------------------
            # this uses a threadpool to compute the jobs
            #--------------------------------------------------------
            
            #print "Using %d threads" % self.threads

            inputs = []

            for molecule, snp_positions in self.snp_positions.iteritems():

                print "Setting up %d SNP jobs for molecule '%s'" % (len(snp_positions),molecule)

                for snp_pos in snp_positions:
                
                        inputs.append((molecule, snp_pos,
                                       self.contigs, self.blast_data, self.genbank_refs,
                                       self.flanking_bases, self.blast_threshold,self.blast_length_threshold))


            if self.threads < 0:

                print "Value of less than zero (%d) in threads, calculating to use 30%% of %d SNPs" % (self.threads,len(inputs))
                self.threads = len(inputs)/3
                
            print "Running %d individual input jobs for processing on %d threads..." % (len(inputs), self.threads)

            pool = ThreadPool(self.threads)
            pool_output = pool.map(compute_snp, inputs)
            pool.close()
            pool.join()
            
            print "Collecting results..."

            # clearing out this memory
            inputs = None

            self.verified_snps = [r for r in pool_output]
             
        else:
            # this is the regular serial mode. Simplest code-wise.

            for molecule, snp_positions in self.snp_positions.iteritems():

                print "Checking %d SNP positions for '%s'" % (len(snp_positions),molecule)

                for snp_pos in snp_positions:
                
                    try:
                    
                        # this method creates the snp and adds it to the global list of snps
                        self.verified_snps.append(SNP( molecule, snp_pos,
                                                       self.contigs, self.blast_data, self.genbank_refs,
                                                       self.flanking_bases, self.blast_threshold, self.blast_length_threshold))
                        
                    except Exception ,e:

                        print "Unable to create a snp entry for molecule '%s' at snp position %d: %s" % (molecule, snp_pos, e)
                        continue

        #- End else clause ------------------------------------------------------------------


        #------------------------------------------------------------------------------------
        # here we count up the snps per gene and the syn/nsyn per gene. 
        # 
        #-------------------------------------------------------------------------------------
        print "\nCounting the snps per gene for %d snps\n" % len(self.verified_snps)

        for s in self.verified_snps:

            if s.total_hits > 0:

                # add a snp for this gene
                self.add_snp_per_gene(s.molecule, s.gene_name)

                # here we add all syns and nsyns per gene by looping through all
                # collected.
                for stat in s.all_stats:

                    self.add_syn_nsyn_per_gene(s.molecule, s.gene_name, stat)
                        
    #-------------------------------------------------------------------------------        

    def write_table(self, output_file):

        # here we write all of the snps and queries to the table
        outfile = open(output_file, 'w')

        print "Writing merged table to '%s'" % output_file

        outfile.write(self.table_header()+"\n")


        for s in self.verified_snps:

            outfile.write(self.snp_to_string(s)+"\n")
            
        outfile.close()
        
    #-------------------------------------------------------------------------------

    def get_snps_per_gene(self, molecule, gene):

        try:

            return str(self.snps_per_gene[molecule][gene])

        except KeyError:

            if gene == "intergenic":
                
                return "intergenic"

            else:

                return "0"
            
    #-------------------------------------------------------------------------------

    def get_snp_div_gene_length(self, snp):

        snps_p_gene = self.get_snps_per_gene(snp.molecule_name, snp.gene_name)

        if snps_p_gene == "intergenic":
            return "intergenic"

        else:

            spg = int(snps_p_gene)

            return str( float(spg)/float(snp.gene_length) )

    #-------------------------------------------------------------------------------

    def get_syn_to_nsyn_ratio(self, snp):

        if snp.gene_name == "intergenic":
            return "intergenic"

        syns = 0
        nsyns = 0
        try:
            syns = self.syn_nsyn_per_gene[snp.molecule_name][snp.gene_name]['SYN']
        except KeyError:
            syns = 0

        try:
            nsyns = self.syn_nsyn_per_gene[snp.molecule_name][snp.gene_name]['NSYN']
        except KeyError:
            nsyns = 0

        if syns == 0:
            return "%d:%d" % (nsyns,syns)
        else:
            return str( float(nsyns)/float(syns) )

    #-------------------------------------------------------------------------------
    # Adds one snp per gene on the molecule to the class-wide dict. One call will add one
    # snp to the count. 
    #-------------------------------------------------------------------------------
    def add_snp_per_gene(self, molecule, gene):
            
        #------------------------------------------------------------
        def increment_gene(gene):

            if self.snps_per_gene[molecule].has_key(gene):

                self.snps_per_gene[molecule][gene] += 1
                
            else:

                self.snps_per_gene[molecule][gene] = 1
        #------------------------------------------------------------

        # if it's intergenic we don't count it here
        # probably doesn't get here but doing it just in case.
        if gene == "intergenic":

            return

        if self.snps_per_gene.has_key(molecule):

                increment_gene(gene)

        else:

            self.snps_per_gene[molecule] = {gene : 1} 


    #-------------------------------------------------------------------------------
    # subtracts one from the gene count. Does nothing if the keys are not present.
    #-------------------------------------------------------------------------------
    def subtract_snp_per_gene(self, molecule, gene):

        #------------------------------------------------------------
        def decrement_gene(gene):

            if self.snps_per_gene[molecule].has_key(gene):

                if self.snps_per_gene[molecule][gene] == 0:
                    return
                else:
                    self.snps_per_gene[molecule][gene] -= 1
                
        #------------------------------------------------------------

        if gene == "intergenic":
            return

        if self.snps_per_gene.has_key(molecule):

                decrement_gene(gene)


    #-------------------------------------------------------------------------------
    # 
    #-------------------------------------------------------------------------------
    def add_syn_nsyn_per_gene(self, molecule, gene, stat):
            
        #------------------------------------------------------------
        def increment_stat(st):

            if self.syn_nsyn_per_gene[molecule][gene].has_key(st):

                self.syn_nsyn_per_gene[molecule][gene][st] += 1
                
            else:

                self.syn_nsyn_per_gene[molecule][gene][st] = 1
        #------------------------------------------------------------

        # if it's intergenic we don't count it here
        # probably doesn't get here but doing it just in case.
        if gene == "intergenic" or stat in ['--', 'No Hit','indel']:

            return

        if self.syn_nsyn_per_gene.has_key(molecule):

            if self.syn_nsyn_per_gene[molecule].has_key(gene):

                increment_stat(stat)

            else:

                if stat == "SYN":
                    self.syn_nsyn_per_gene[molecule][gene] = {'SYN' : 1, 'NSYN': 0}
                elif stat == "NSYN":
                    self.syn_nsyn_per_gene[molecule][gene] =  {'SYN' : 0, 'NSYN': 1}
                else:
                    raise Exception("Invalid stat value to add to the ds/dn: %s" % stat)


        else:
            if stat == "SYN":
                self.syn_nsyn_per_gene[molecule] = {gene : {'SYN' : 1, 'NSYN': 0}}
            elif stat == "NSYN":
                self.syn_nsyn_per_gene[molecule] = {gene : {'SYN' : 0, 'NSYN': 1}}
            else:
                raise Exception("Invalid stat value to add to the ds/dn: %s" % stat)
                
            
    #-------------------------------------------------------------------------------------------

    def parse_genbank_list(self):
        """
        This method parses a file that contains a list of genbank files, one per line. It opens each and loads
        it into a list.
        """
        
        num_gens = 0
        genbank_files = []

        print "Loading Genbank files"
        
        if not self.genbank_list_file is None:

            g_files = self.parse_list_file(self.genbank_list_file)
            
            print "\tGot %d Genbank files from list file '%s'" % (len(g_files), self.genbank_list_file)

            genbank_files += g_files
            

        if not self.genbank_list is None:

            print "\tGot %d Genbank files from the command line" % len(self.genbank_list)

            genbank_files += self.genbank_list
           

        # now we can loop through all of the files and load them
        for gbf in genbank_files:

            try: 
                gbdata = self._parse_genbank_file_list(gbf)
            except Exception, e:
                print "Can't parse genbank file %s: %s" % (gbf, e)
                continue
            
            self.genbank_refs.append(gbdata)

                
    #-------------------------------------------------------------------------------
    def _parse_genbank_file(self,genbank_file):
        
        genbank_recs = []
        num_gb_recs = 0

        genbank_input_handle = open(genbank_file, "rU")
        for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

            genbank_recs[gbrec.id] = gbrec
            num_gb_recs += 1

        genbank_input_handle.close()

        print "Read in %i genbank records from file %s" % (num_gb_recs, genbank_file)

        return genbank_recs

    #-------------------------------------------------------------------------------

    def _parse_genbank_file_dict(self,genbank_file):
        
        genbank_recs = {}
        num_gb_recs = 0

        genbank_input_handle = open(genbank_file, "rU")
        for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

            genbank_recs[gbrec.name] = gbrec
            num_gb_recs += 1

        genbank_input_handle.close()
        genbank_recs
        print "Read in %i genbank records from file %s" % (num_gb_recs, genbank_file)

        return genbank_recs

    #-------------------------------------------------------------------------------

    def _parse_genbank_file_list(self,genbank_file):
        
        genbank_recs = []

        genbank_input_handle = open(genbank_file, "rU")
        for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

            genbank_recs.append(gbrec)

        genbank_input_handle.close()

        print "Read in %i genbank records from file %s" % (len(genbank_recs), genbank_file)

        return genbank_recs

    #-------------------------------------------------------------------------------


    def parse_snp_positions(self):


        num_snp_filess = 0
        snp_panel_files = []

        print "Loading SNP panel files"
        
        if not self.snp_list_file is None:

            s_files = self.parse_list_file(self.snp_list_file)
            
            print "\tGot %d SNP panel files from list file '%s'" % (len(s_files), self.snp_list_file)

            snp_panel_files += s_files
            

        if not self.snp_list is None:

            print "\tGot %d SNP panel files from the command line" % len(self.snp_list)

            snp_panel_files += self.snp_list


        
        for s in snp_panel_files:

            if os.path.isfile( os.path.abspath(s) ):

                self._parse_snp_positions_dict(os.path.abspath(s))

            else:

                print "Can't parse '%s', it's not a valid file" % s
            
        
        for key, value in self.snp_positions.items():

                self.snp_positions[key].sort()


        print "\n"

    #-------------------------------------------------------------------------------
    # This does the same as above but returns a bucketed dict 
    def _parse_snp_positions_dict(self, snp_file):

        num_snps = 0

        with open(snp_file, "rU") as input_handle:

            for line in input_handle:

                if not line.isspace():

                    snp_txt = line.split()

                    if len(snp_txt) != 2:
                        print "Error parsing entry \"%s\", not in form 'locus position'" % line
                        continue

                    else:

                        if self.snp_positions.get(snp_txt[0]) is None:

                            try:
                                self.snp_positions[snp_txt[0]] = [int(snp_txt[1])]
                            except Exception,e:
                                raise Exception("Error with snp: %s" % e)

                        else:

                            # if this refpos is in the positions for this locus, when we ignore it,
                            # otherwise It needs to be added.
                            try:
                                ref_pos_value = int(snp_txt[1])
                            except Exception, e:
                                print "Can't add '%s' to locus %s: %s" % (snp_txt[1], snp_txt[0], e)
                                continue
                                    
                            if not ref_pos_value in self.snp_positions[snp_txt[0]]:

                               
                                    self.snp_positions[snp_txt[0]].append(ref_pos_value)

                            else:

                                print "SNP position %s is already in %s, skipping" % (snp_txt[1], snp_txt[0])
                    
                    num_snps += 1

            
        
        print "Parsed %d snp positions from file '%s'" % (num_snps,snp_file)

    #-------------------------------------------------------------------------------
    # Using some example code I found from this tutorial
    #
    #   http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec100
    #
    # Might need to use this:
    #  http://biopython.org/DIST/docs/api/Bio.SearchIO.BlastIO-module.html
    #
    #-------------------------------------------------------------------------------
    def parse_blast_data(self):

        num_blast_files = 0

        blast_files = []

        print "Loading BLAST result files"
        
        if not self.blast_list_file is None:

            b_files = self.parse_list_file(self.blast_list_file)
            
            print "\tGot %d BLAST result files from list file '%s'" % (len(b_files), self.blast_list_file)

            blast_files += b_files
            

        if not self.blast_list is None:

            print "\tGot %d BLAST result files from the command line" % len(self.blast_list)

            blast_files += self.blast_list


        print "Parsing %d BLAST result files" % len(blast_files)


        for blast_file in blast_files:
            
            filename = blast_file.rstrip()
                
            if os.path.isfile( os.path.abspath(filename) ):

                read_successful = True
                try:
                    blast_rec = self._parse_blast_dict(filename) # loads this as a list of lists
                except Exception, e:

                    print "Can't parse '%s': %s" % (filename, e)
                    read_successful = False
                    continue

                # we only out it and report success if the parse_blast method
                # is successful.
                if read_successful:

                    print "Parsed BLAST result file: %s" % filename
                    
                    num_blast_files += 1

            else:

                print "Can't parse BLAST result file '%s', not a valid file" % filename


        print "Parsed %d BLAST result files\n" % num_blast_files

        if self.blast_threshold > 0.0:

            num_hits_filtered = self.num_blast_hits_loaded - self.num_filtered_blast_hits
            
            print "Filtered out %d blast hits of %d total with %f threshold\n" % (num_hits_filtered, self.num_blast_hits_loaded, self.blast_threshold)

    #-------------------------------------------------------------------------------

    def parse_list_file(self, list_file):
        """
        Just parses a list file with a list of files one per line.
        """

        files = []
        
        with open(list_file, "rU") as input_handle:

            for f in input_handle:

                line = f.rstrip()

                if not line == "":
                    files.append(line)


        return files

                
    #-------------------------------------------------------------------------------

    def parse_query_genomes(self):

        num_queries = 0
        num_files = 0

        query_contig_files = []


        print "Loading query contig FASTA files"
        
        if not self.query_contig_list_file is None:

            qc_files = self.parse_list_file(self.query_contig_list_file)
            
            print "\tGot %d query contig genome files from list file '%s'" % (len(qc_files), self.query_contig_list_file)

            query_contig_files += qc_files
            

        if not self.query_contig_list is None:

            print "\tGot %d query contig genome files from the command line" % len(self.query_contig_list)

            query_contig_files += self.query_contig_list


        print "Parsing %d query contig genome files" % len(query_contig_files)

        # loop through each file and load them
        for qc_file in query_contig_files:

            self.contigs.append( ContigFile(qc_file) )

        print "Parsed %d query contig genome files\n" % len(self.contigs)


#-------------------------------------------------------------------------------

    def _filter_query_result(self, qr):
        """
        Returns a tuple with the query id and filtered hits.
        """
        bitscore_filter = lambda hsp: hsp.bitscore >= self.blast_threshold 

        query_id = qr.id
        
        filtered_hits = []

        # collect all hits that pass the bitscore filter
        for h in qr.hits:

            fh = h.filter(bitscore_filter)

            if not fh is None:
                filtered_hits.append(h)

        return (query_id, filtered_hits)
    
#-------------------------------------------------------------------------------

    def _parse_blast_dict(self, blast_file):
        """
        full example for parsing blast files here:
            http://bugzilla.open-bio.org/attachment.cgi?id=293&action=view
        """

        for qresult in SearchIO.parse(blast_file, 'blast-text'):

            qid = None
            hits = []

            # getting number of individual blast hits and keeping track
            self.num_blast_hits_loaded += sum([len(h) for h in qresult.hits])


            if self.blast_threshold > 0.0:

                (qid, hits) = self._filter_query_result(qresult)

                # adds up the number of filtered hits
                self.num_filtered_blast_hits += sum([len(h) for h in hits])
                

            qid = qresult.id
            hits = qresult.hits

            #------------------------------------------------------
            # If we're doing multiprocessing we need to make the QueryResult
            # pickleable for the forking. It may also cut down on overall
            # memory size. Also since I don't want to bother making a new QueryResult
            # this is easier to keep around.
            #------------------------------------------------------
            qresult = PickleableQueryResult(qid, hits)

            if self.blast_data.has_key(qid):
                
                self.blast_data[qid].append(qresult)

            else:

                self.blast_data[qid] = [qresult]


#-------------------------------------------------------------------------------
def __main__():


    #Parse Command Line
    parser = argparse.ArgumentParser()

    
    blast_group = parser.add_argument_group(title="BLAST result files", description="Raw text output from a blastn run")
    blast_group.add_argument('-w', '--blast-list-file', help="a file with a list to blast files to parse, with one per line", metavar="bfiles.list")
    blast_group.add_argument('-b', '--blast', nargs='*', help="a list of raw blast files separated by spaces", metavar="bfile.blastn")

    genbank_group = parser.add_argument_group(title="Genbank files", description="Genbank files")
    genbank_group.add_argument('-r', '--genbank-list-file', help="A file that contains a list of genbank files with one per line", metavar="gbfiles.list")
    genbank_group.add_argument('-g', '--genbank', nargs='*', help="list of genbank annotation files separated by spaces", metavar="gbfile.gb")

    snp_panel_group = parser.add_argument_group(title="SNP panel files", description="Contains a molecule name that maps to the Genbank annotation along with a SNP position, one per line")
    snp_panel_group.add_argument('-y', '--snp-panel-list-file', help="A file that contains a list of snp panel files with one per line", metavar="snp_panel.list")
    snp_panel_group.add_argument('-s', '--snp-panel', nargs='*', help="list of snp position files separated by spaces", metavar="snp_panel.tabular")

    query_contig_group = parser.add_argument_group(title="Query Contig FASTA files", description="FASTA files that contain extracted regions from a genome. Query Contig files can be listed with a metaname appended via a colon 'query_file:meta_name'. The meta_name will be present in the table column header.")
    query_contig_group.add_argument('-z', '--query-list-file', help="File with list of query genomes with one per line", metavar="qfiles.list")
    query_contig_group.add_argument('-q', '--query', nargs='*', help="list of query files separated by spaces.", metavar="qfile.fasta")

    
    parser.add_argument('-f','--flanking-bases',  type=int, required=True, help='Number of bases on *each side* of the SNP positions in the reference genome.' )
    parser.add_argument('-o', '--out', help="Output file for the merged table", default="merged_table.txt", metavar="merged_table.txt")
    parser.add_argument('-t', '--blast-threshold', type=str, help="Takes a percentage for a blast score as an int or decimal: ex 80%% can be 80 or 0.8")
    parser.add_argument('-l', '--blast-length-threshold', type=int, default=0, help="The minimum length of a blast hit to count as a SNP")
    parser.add_argument('-a','--amino-table',  type=int, default=1, help="The amino table to use when translating codons. Default is 1." )
    parser.add_argument('-v', '--version', action='version', version=__version__)
    

    procesing_group = parser.add_mutually_exclusive_group()
    procesing_group.add_argument('-m', '--multiproc', type=int, help="Number of processes to use")
    procesing_group.add_argument('-n', '--threads', type=int, help="Number of threads to use")


    args = parser.parse_args()


    if args.blast_list_file is None and args.blast is None:

        parser.print_usage()
        print "You need at least a blast list file containing a listing with one blast file per line with the -b option, or to have at least one blast file listed with the -w option."
        print "\n"
        
        sys.exit()

    if args.genbank_list_file is None and args.genbank is None:

        parser.print_usage()
        print "You need at least a genbank list file containing a listing with one genbank file per line with the -r option, or to have at least one genbank file listed with the -x option."
        print "\n"
        sys.exit()

    if args.query_list_file is None and args.query is None:

        parser.print_usage()
        print "You need at least a query list file containing a listing with one query contig fasta file per line with the -q option, or to have at least one query contig fasta file listed with the -z option."
        print "\n"
        sys.exit()

    if args.snp_panel_list_file is None and args.snp_panel is None:

        parser.print_usage()
        print "You need at least a SNP panel list file containing a listing with one SNP panel file per line with the -s option, or to have at least one SNP panel file listed with the -y option."
        print "\n"
        sys.exit()



    start_time = timeit.default_timer()
    
    snp_verify = SNPVerify(blast_list_file=args.blast_list_file,
                           genbank_list_file=args.genbank_list_file,
                           snp_list_file=args.snp_panel_list_file,
                           query_contig_list_file=args.query_list_file,
                           flanking_bases=args.flanking_bases,
                           amino_table=args.amino_table,
                           blast_list=args.blast,
                           genbank_list=args.genbank,
                           snp_list=args.snp_panel,
                           query_contig_list=args.query,
                           blast_threshold=args.blast_threshold,
                           blast_length_threshold=args.blast_length_threshold,
                           multiproc=args.multiproc,
                           threads=args.threads,
                           )

    # write out the table 
    snp_verify.write_table(args.out)

    elapsed_time = timeit.default_timer() - start_time

    mins, secs = divmod(elapsed_time, 60)
    hours, mins = divmod(mins, 60)
    days, hours = divmod(hours, 24)

    print "Execution time was: %dD:%02dH:%02dM:%02dS" % (days, hours, mins, secs)

    print "Done!"

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
