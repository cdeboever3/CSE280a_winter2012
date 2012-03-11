# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os, re
import numpy as np

### helper functions ###

p_match = re.compile('\.|,') # regular expression that matches characters that indicate a match in the pileup file
p_alt = re.compile('[ATCG]',flags=re.IGNORECASE) # regex that matches alternate bases in the pileup file

def parse_plineL(L):
    """
    Parse pileup line into output format
    Input: Pileup line that has been split by \t into a list
    Output: Separated by \t: chrom, position, ref base, alt base, number of reads supporting reference, number of reads supporting alternate
    """
    refL = p_match.findall(L[4])
    altL = p_alt.findall(L[4])
    ref = len(refL)
    alt = len(altL)
    return '\t'.join([lineL[0],(str(int(lineL[1])-1)),lineL[2],altL[0],str(ref),str(alt)]) # subtract one from lineL[1] to make 0-based (pileup is 1-based) coordinate
    return '{0}\t{1}'.format(ref,alt)

def filter_site(L,cc):
    """
    Decide whether pileup line should be skipped according to various filters
    Input: Pileup line that has been split by \t into a list and coverage cutoff
    Output: Boolean indicating whether pileup line should be skipped
    """
    skip = False
    if lineL[2] == 'N': # skip if reference base was an N
        skip = True
    if '+' in lineL[4] or '-' in lineL[4] or '^' in lineL[4] or 'N' in lineL[4]: # skip if there is an indel or an uncalled base (indicative of bad quality site)
        skip = True
    if int(3) < cc: # skip if the read coverage of the base is below the cutoff
        skip = True
    if len(set(p_alt.findall(L[4]))) > 1: # skip if there is more than one alternate base at the site (indicative of bad quality)
        skip = True 
    if len(p_alt.findall(L[4])) == 0: # skip if there is no alternate allele present
        skip = True
    return skip

if __name__ == '__main__':

    ### magic variables ###
    # these variables can be set at the command line as well
    cov_cutoff = 5 # coverage cutoff for each site

    ### gather command line arguments ###
    parser = argparse.ArgumentParser(description='This script takes takes a pileup file and parses the pileup file to make the tsv file needed for mixed_variant_calling.py.')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('-c', metavar='cov_cutoff', default=cov_cutoff, type=int, help='Coverage cutoff for each site. If the coverage for a site is below this threshold, the site will not be included.')
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    inF         = args.infile
    outF        = args.outfile
    cov_cutoff  = args.c
    debug       = args.d

    print >>outF, '#chromosome\tposition (0-based)\tref base\talt base\tref coverage\talt coverage'

    line = inF.readline().strip()

    while line != '':
        lineL = line.split('\t')
        skip = filter_site(lineL,cov_cutoff)
        if not skip:
            print >>outF, parse_plineL(lineL)
        line = inF.readline().strip()
