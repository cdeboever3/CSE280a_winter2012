# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os
import numpy as np
from bisect import bisect_left 
from scipy.stats import binom

### helper functions ###

def find_lt(a,x):
    """
    Find rightmost value less than x in list a
    Input: list a and value x
    Output: rightmost value less than x in a
    """
    i = bisect_left(a,x)
    if i:
        return a[i-1]
    raise ValueError

def find_ge(a,x):
    """
    Find leftmost item greater than or equal to x in list a
    Input: list a and value x
    Output: leftmost value less than or equal to x in a
    """
    i = bisect_left(a,x)
    if i != len(a):
        return a[i]
    raise ValueError

def get_altL(fn):
    """
    Make a list of alternate allele frequencies and number of reads
    Input: tsv file with reference freq in first column and alterate freq in second column
    Output: a list of tuples with number of reads and alternate allele frequency
    """
    f = open(fn,'r')
    linesL = [ x.strip().split('\t') for x in f.readlines() ]
    f.close()
    for i in range(len(linesL)):
        if linesL[i][0] == '0':
            linesL[i][0] = '1'
    return zip([ int(x[0])+int(x[1]) for x in linesL ], [ float(x[1])/(float(x[0])+float(x[1])) for x in linesL ]) # each tuple is [freq,num_reads]

# def generate_cancer_possible_freqL(pL,sL,er): 
# I want to make a function which generates the likely frequencies seen in a cancer sample. This would exclude double-hit mutations (i.e. a single site gains somatic mutations on both chromosomes). This simplifications can only be made in the diploid case, however, because ploidy-variable populations might be weird...

def generate_possible_freqL(pL,sL,er): 
    """
    Generate list of possible allele frequencies
    Input: ploidy list, frequency (of each subpopulation) list, and sequencing error rate
    Output: list of possible allele frequences
    """
    h = sum(pL) # number of different haplotypes
    L = [ bin(x)[2:] for x in range(0,2**h) ]
    M = [ '0'*(len(L[-1])-len(x))+x for x in L ]
    p_freqL = []
    for i in range(len(pL)):
        p_freqL += [sL[i]/pL[i]]*pL[i]
    p_freqA = np.array(p_freqL)
    sA = np.array(sL)
    aL = []
    for g in M:
        aL.append(sum(np.array([ int(x) for x in list(g) ])*p_freqL))
    freqL = sorted(list(set(aL))) # only look at distinct frequencies
    freqL[0] = er # we don't want 0 for binomial pmf or we'll end up with 0 probabilities
    freqL[-1] = 1-er # we don't want 1 for binomial pmf or we'll end up with 0 probabilities
    return freqL
    
def grid_search_parameters(step):
    """
    Make a list of parameters to try
    Input: step size
    Output: subpopulation frequencies to try
    """
    f1 = list(np.arange(0,1+step,step))
    f2 = list(np.arange(0,1+step,step))
    f2.reverse()
    return zip(f1,f2)

if __name__ == '__main__':

    ### magic variables ###
    # these variables can be set at the command line as well
    ploidyL = [2,2] # the entries in this list are the expected ploidy of each subpopulation. Default is two diploid subpopulations
    error_rate = 0.001 # sequencing error rate

    ### gather command line arguments ###
    parser = argparse.ArgumentParser(description='This script determines the relative frequencies of different populations and estimates the genotypes.')
    parser.add_argument('infile', help='Input tsv file. First column is number of reads supporting reference allele and second column is number of reads supporting alternate allele.')
    parser.add_argument('-pL', default=ploidyL, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-er', default=error_rate, type=float, help='Sequencing error rate. For instance, 0.01 means that 1/100 base calls will be incorrect. Default: {0}'.format(error_rate))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    inN         = args.infile
    ploidyL     = args.pL
    error_rate  = args.er
    debug       = args.d

    inN = os.path.realpath(inN) # get the input file path

    if len(ploidyL) > 2:
        print >>sys.stderr, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    altL = get_altL(inN) # a list of number of reads and alternate allele frequencies

    ### grid search ###

    parL = grid_search_parameters(0.01)
    best_par = []
    best_ll = float("-inf")

    for par in parL:
        exp_freqL = generate_possible_freqL(ploidyL,par,error_rate)

        ll = 0 # log-likelihood

        for alt in altL:
            try:
                i = find_lt(exp_freqL,alt[1]) # Find rightmost value less than x
            except ValueError:
                i = float("-inf")
            try:
                j = find_ge(exp_freqL,alt[1]) # Find leftmost item greater than or equal to x
            except ValueError:
                j = float("inf")
            if alt[1]-i < j-alt[1]:
                exp_freq = i
            else:
                exp_freq = j
            
            ll += np.log(binom.pmf(round(alt[0]*alt[1]),alt[0],exp_freq))
        
        if ll > best_ll:
            best_ll = ll
            best_par = par

    print >>sys.stdout, 'log-likelihood: {0}\npopulation frequencies: {1}'.format(best_ll,' '.join([ str(x) for x in best_par ]))
