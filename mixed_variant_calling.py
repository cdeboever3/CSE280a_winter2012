# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os

### helper functions ###

def get_alt_freqL(fn):
    """
    Make a list of alternate allele frequencies
    Input: tsv file with reference freq in first column and alterate freq in second column
    Output: a list of alternate allele frequencies
    """
    f = open(fn,'r')
    linesL = [ x.strip().split('\t') for x in f.readlines() ]
    f.close()
    return [ float(x[1])/float(x[0]) for x in linesL ] 

def generate_possible_freqL(pL,aL):
    """
    Generate list of possible allele frequencies
    Input: ploidy list and abundance (of each subpopulation) list
    Output: list of possible allele frequences
    """
    
def grid_search_parameters(step):
    """
    Make a list of parameters to try
    Input: step size
    Output: abundances to try
    """
    f1 = list(np.arange(0,1+step,step))
    f2 = list(np.arange(0,1+step,step))
    f2.reverse()
    return zip(f1,f2)

if __name__ == '__main__':

    ploidyL = [2,2] # the entries in this list are the expected ploidy of each subpopulation. Default is two diploid subpopulations

    parser = argparse.ArgumentParser(description='This script determines the relative abundance of different populations and estimates the genotypes.')
    parser.add_argument('infile', help='Input tsv file. First column is number of reads supporting reference allele and second column is number of reads supporting alternate allele.')
    parser.add_argument('-pL', default=ploidyL, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    inN     = args.infile
    ploidyL = args.pL
    debug   = args.d

    inN = os.path.realpath(inN) # get the input file path

    alt_freqL = get_alt_freqL(inN) # a list of alternate allele frequencies

    if len(ploidyL) > 2:
        print >>sys.err, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    parL = grid_search_parameters(0.01)
