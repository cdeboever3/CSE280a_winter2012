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
    if lineL[0][0] == '#':
        lineL = lineL[1:]
    for i in range(len(linesL)):
        if linesL[i][4] == '0': # if the number of reads supporting alternate is 0, we'll switch to 1 so avoid numeric issues
            linesL[i][4] = '1'
    return zip([ int(x[4])+int(x[5]) for x in linesL ], [ float(x[5])/(float(x[4])+float(x[5])) for x in linesL ]) # each tuple is [freq,num_reads]

# def generate_cancer_possible_freqL(pL,sL,er): 
# I want to make a function which generates the likely frequencies seen in a cancer sample. This would exclude double-hit mutations (i.e. a single site gains somatic mutations on both chromosomes). This simplifications can only be made in the diploid case, however, because ploidy-variable populations might be weird...

def generate_possible_freqL(pL,sL,er): 
    """
    Generate list of possible allele frequencies
    Input: ploidy list, frequency (of each subpopulation) list, and sequencing error rate
    Output: list of possible allele frequences
    """
    h = sum(pL) # number of different haplotypes
    L = [ bin(x)[2:] for x in range(1,2**h-1) ] # range from 1 to 2**h-1 because we don't want 0% or 100% allele freq
    M = [ '0'*(len(L[-1])-len(x))+x for x in L ]
    p_freqL = []
    for i in range(len(pL)):
        p_freqL += [sL[i]/pL[i]]*pL[i]
    p_freqA = np.array(p_freqL)
    sA = np.array(sL)
    aL = []
    for g in M:
        aL.append(sum(np.array([ int(x) for x in list(g) ])*p_freqL))
    return sorted(list(set(aL+[er,1-er]))) 

def freq_to_genotype(pL,sL,er): 
    """
    Creates dict of expected alternate allele frequencies and consistent genotypes
    Input: ploidy list, frequency (of each subpopulation) list, and sequencing error rate
    Output: dict of expected alternate allele frequencies and consistent genotypes. Genotypes represented as binary strings in the order of the ploidy list
    """
    h = sum(pL) # number of different haplotypes
    L = [ bin(x)[2:] for x in range(1,2**h-1) ] # range from 1 to 2**h-1 because we don't want 0% or 100% allele freq
    M = [ '0'*(len(L[-1])-len(x))+x for x in L ]
    p_freqL = []
    for i in range(len(pL)):
        p_freqL += [sL[i]/pL[i]]*pL[i]
    p_freqA = np.array(p_freqL)
    sA = np.array(sL)
    aD = {} # dict where each key is an expected alternate allele frequency and each value is a list of genotypes consistent with this alternate allele frequency
    for g in M:
        alt_freq = sum(np.array([ int(x) for x in list(g) ])*p_freqL)
        if aD.has_key(g):
            aD[alt_freq].append(g)
        else:
            aD[alt_freq] = [g]
    return aD
    
def grid_search_parameters(step):
    """
    Make a list of parameters to try
    Input: step size
    Output: subpopulation frequencies to try
    """
    f1 = list(np.arange(step,1,step))
    f2 = list(np.arange(step,1,step))
    f2.reverse()
    return zip(f1,f2)

def estimate_genotype(alt_freq,exp_freqL):
    """
    Maximum likelihood estimator of alt_freq given possibilities in exp_freqL
    Input: observed alternate frequency and list of expected alternate frequencies
    Output: ML estimator of true alternate allele frequency
    """
    try:
        i = find_lt(exp_freqL,alt_freq) # Find rightmost value less than x
    except ValueError:
        i = float("-inf")
    try:
        j = find_ge(exp_freqL,alt_freq) # Find leftmost item greater than or equal to x
    except ValueError:
        j = float("inf")
    if alt_freq-i < j-alt_freq:
        return i
    else:
        return j

def main():
    ### magic variables ###
    # these variables can be set at the command line as well
    ploidyL = [2,2] # the entries in this list are the expected ploidy of each subpopulation. Default is two diploid subpopulations
    error_rate = 0.001 # sequencing error rate

    ### gather command line arguments ###
    parser = argparse.ArgumentParser(description='This script determines the relative frequencies of different populations and estimates the genotypes.')
    parser.add_argument('infile', help='Input tsv file. Columns should be: chrom, position, ref base, alt base, number of reads supporting reference, number of reads supporting alternate.')
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

    ### find population frequencies ###

    parL = grid_search_parameters(0.01) # grid search
    best_par = []
    best_ll = float("-inf")

    for par in parL:
        exp_freqL = generate_possible_freqL(ploidyL,par,error_rate)
        ll = 0 # log-likelihood

        for alt in altL:
            exp_freq = estimate_genotype(alt[1],exp_freqL)
            ll += np.log(binom.pmf(round(alt[0]*alt[1]),alt[0],exp_freq)) 
            # round(alt[0]*alt[1]) is the number of reads we saw supporting alternate allele (i.e. the number of successes under the binomial test)
            # alt[0] is the total number of reads covering this site (i.e. the number of attempts in our binomial test)
            # exp_freq is our probability of success (i.e. observing a read supporting alternate) from our ML estimation (see estimate_genotype)
        
        if ll > best_ll:
            best_ll = ll
            best_par = par

    ### determine genotypes ###
    # use best population frequency parameters and walk through sites, assign genotypes, p-values or scores maybe?

    print >>sys.stdout, 'log-likelihood\t{0}\npopulation frequencies\t{1}'.format(best_ll,'\t'.join([ str(x) for x in best_par ]))

if __name__ == '__main__':
    main()
