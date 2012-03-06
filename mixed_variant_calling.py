# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os
import numpy as np

### helper functions ###

def find_lt(a,x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def get_altL(fn):
    """
    Make a list of alternate allele frequencies
    Input: tsv file with reference freq in first column and alterate freq in second column
    Output: a list of tuples with number of reads alternate allele frequency
    """
    f = open(fn,'r')
    linesL = [ x.strip().split('\t') for x in f.readlines() ]
    f.close()
    for i in range(len(linesL)):
        if linesL[i][0] == '0':
            linesL[i][0] = '1'
    return zip([int(x[0])+int(x[1]) for x in linesL ], [ float(x[1])/float(x[0]) for x in linesL ])

def generate_possible_freqL(pL,sL):
    """
    Generate list of possible allele frequencies
    Input: ploidy list and frequency (of each subpopulation) list
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
    return sorted(list(set(aL))) # only look at distinct frequencies

    
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

    ploidyL = [2,2] # the entries in this list are the expected ploidy of each subpopulation. Default is two diploid subpopulations

    parser = argparse.ArgumentParser(description='This script determines the relative frequencies of different populations and estimates the genotypes.')
    parser.add_argument('infile', help='Input tsv file. First column is number of reads supporting reference allele and second column is number of reads supporting alternate allele.')
    parser.add_argument('-pL', default=ploidyL, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    inN     = args.infile
    ploidyL = args.pL
    debug   = args.d

    inN = os.path.realpath(inN) # get the input file path

    altL = get_altL(inN) # a list of number of reads and alternate allele frequencies

    if len(ploidyL) > 2:
        print >>sys.stderr, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    ### grid search ###

    parL = grid_search_parameters(0.01)
    best_parL = []
    best_ll = float("inf")

    for par in parL:
        exp_freqL = generate_possible_freqL(ploidyL,parL)

        ll = 0 # log-likelihood

        for alt in altL:
            i = find_lt(exp_freqL,alt[1])
            j = find_ge(exp_freqL,alt[1])
            if alt[1]-exp_freqL[i] < exp_freqL[j]-alt[1]:
                exp_freq = exp_freqL[i]
            else:
                exp_freq = exp_freqL[j]

            ll += np.log(binom.pmf(alt[0]*alt[1],alt[0],exp_freq))
        
        if ll < best_ll:
            best_ll = ll
            best_parL = parL

    print >>sys.stdout, 'log-likelihood: {0}\npopulation frequencies: {1}'.format(ll,' '.join(best_parL))

