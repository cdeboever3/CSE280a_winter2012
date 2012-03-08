# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os, subprocess
import numpy as np

### helper functions ###

def parse_plineL(L):
    """
    Parse pileup line into output format
    Input: Pileup line that has been split by \t into a list
    Output: Number of reads supporting reference, number of reads supporting alternate, separated by \t
    """
    ref = regex
    alt = regex

if __name__ == '__main__':

    ### magic variables ###
    # these variables can be set at the command line as well
    cov_cutoff = 10 # coverage cutoff for each site

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

    line = f.readline().strip()

    while line != '':
        lineL = line.split('\t')
        skip = False
        if lineL[2] == 'N':
            skip = True
        if '+' in lineL[4] or '-' in lineL[4] or '^' in lineL[4]:
            skip = True
        if int(3) < cov_cutoff:
            skip = True

        if not skip:
            print >>outF, parse_plineL(lineL)

        line = f.readline().strip()
