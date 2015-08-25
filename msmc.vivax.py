# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:38:37 2015

@author: stsmall
"""

import sys
from collections import defaultdict

#parser = argparse.ArgumentParser(description="This is for converting Fasta to MSMC")
#parser.add_argument('-i', action="append", dest="fasta")

#if default=False, just call to make True; dest links it to name access it via args.dest, e.g., args.fasta

#%% function fastadict
def fasta_dict(fasta):
    """converts fasta to dictionary"""
    msmc_dict = defaultdict(list)#vivax03:AAAAAT,AAAAAT,AAAAAT,AAAAT
    for i in fasta:  #through fasta list as -i
        with open("%s" %i, 'r') as fa1:
            chrom1 = fa1.next()[1:].rstrip('\n')           
            seq = ""
            for line in fa1:
                if '>' not in line:
                    seq += line.rstrip('\n') 
                else:
                    msmc_dict[chrom1].append(seq)                    
                    chrom1 = line[1:].rstrip('\n')
                    seq = ""
            if seq != "":
                msmc_dict[chrom1].append(seq)
    return msmc_dict

#%% function findsnps
def find_snps(msmc_dict):
    """records differences in fasta via dictionary and prints to MSMC inputfile"""
    for chrom in msmc_dict.keys():
        snp_positions = []        
        with open("%s.msmc.in" %chrom, "w") as outfile:
            for ind in msmc_dict[chrom][1:]:                            
                snp_positions += ([i for i in xrange(len(msmc_dict[chrom][0])) if ind[i] != msmc_dict[chrom][0][i]])
            snp_coordinates = sorted(set(snp_positions))
            ncount = []
            phased_seq = ''
            het = 0
            for pos in snp_coordinates:
                for ind in msmc_dict[chrom]:          
                    phased_seq += ind[pos]
                    ncount.append(ind[het:pos].count("N"))
                if "N" in phased_seq:                                   
                    phased_seq = ''
                else:                  
                    homo_run = int(pos)-int(het) - max(ncount)
                    outfile.write("%s\t%s\t%s\t%s\n" %(chrom, pos, homo_run, phased_seq))
                    #print("%s\t%s\t%s\t%s\n" %(chrom, pos, homo_run, phased_seq))                    
                    phased_seq = ''
                    het = pos
                    ncount = []
#%% main call
if __name__ == "__main__":
    find_snps(fasta_dict(sys.argv[1:]))
