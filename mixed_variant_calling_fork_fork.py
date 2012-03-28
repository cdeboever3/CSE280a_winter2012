# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os, re, traceback
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
    if linesL[0][0][0] == '#':
        linesL = linesL[1:]
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
        try:
            p_freqL += [sL[i]/pL[i]]*pL[i]
        except:
            print >>sys.stderr, 'wacka'
            print >>sys.stderr, sL, pL, i
            sys.stderr.flush()
            0 / asdf
    p_freqA = np.array(p_freqL)
    sA = np.array(sL)
    aD = {} # dict where each key is an expected alternate allele frequency and each value is a list of genotypes consistent with this alternate allele frequency
    for g in M:
        alt_freq = sum(np.array([ int(x) for x in list(g) ])*p_freqL)
        if aD.has_key(alt_freq):
            aD[alt_freq].append(g)
        else:
            aD[alt_freq] = [g]
    aD[er] = ['0'*(len(L[-1])-1) + bin(0)[2:]] # add genotype for 0% alternate allele freq
    aD[1-er] = [bin(2**h-1)[2:]] # add genotype for 100% alternate allele freq
    return aD

def collapse_genotypes(pL,gL):
    """
    Reduces a list of genotypes to distinct genotypes given ploidy
    Input: ploidy list pL and list of genotypes gL where each genotype is a binary string ordered according to ploidy list
    Output: genotype list with non-redundant genotypes
    """
    if len(gL) < 2:
        return gL
    else:
        uniqueL = [] # list of unique genotypes relative to ploidy
        for g in gL:
            s = ''
            for i in xrange(len(pL)):
                s += ''.join(sorted(g[0:pL[i]]))
                g = g[pL[i]:]
            if s not in uniqueL:
                uniqueL.append(s)
        return uniqueL
    
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
    cov_cutoff = 4 # coverage cutoff for variant sites

    ### gather command line arguments ###
    parser = argparse.ArgumentParser(description='This script determines the relative frequencies of different populations and estimates the genotypes.')
    parser.add_argument('infile', help='Input tsv file. Columns should be: chrom, position, ref base, alt base, number of reads supporting reference, number of reads supporting alternate.')
    parser.add_argument('-o', nargs='?', type=argparse.FileType('w'),default=sys.stdout, help='Output file. Default: standard out')
    parser.add_argument('-pL', default=ploidyL, type=int, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-er', default=error_rate, type=float, help='Sequencing error rate. For instance, 0.01 means that 1/100 base calls will be incorrect. Default: {0}'.format(error_rate))
    parser.add_argument('-cc', default=cov_cutoff, type=int, help='Coverage cutoff. If the coverage of either the alternate or reference allele is less than or equal to this value, the site will not be considered as a variant site. Default: {0}'.format(cov_cutoff))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    inN         = args.infile
    outF        = args.o
    ploidyL     = args.pL
    error_rate  = args.er
    debug       = args.d
    cov_cutoff  = args.cc

    inN = os.path.realpath(inN) # get the input file path

    if len(ploidyL) > 2:
        print >>sys.stderr, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    ### find population frequencies ###

    parL = grid_search_parameters(0.01) # grid search

    germ_mu = 0.001
    soma_mu = 0.00001

    at_prob = 1/2.  # Marginal probability of normal being AT
    tt_prob = 1/2.  # Marginal probability of normal being TT

    # List of possible genotype pairs (normal, cancer)
    genotypes = [('AA', 'AA'), ('AA','AT'), ('AT','AT'), ('TT','TT')]

    genL = [(0, 0, (1 - germ_mu)*(1 - soma_mu)),                  # ('AA','AA')
            (0, 0.5, (1 - germ_mu) * soma_mu),                    # ('AA','AT')
            (0.5, 0.5, germ_mu * at_prob * (1 - soma_mu)) ,       # ('AT','AT')
            (1, 1, germ_mu * tt_prob * (1 - soma_mu)) ,           # ('TT','TT')
             ]

    np.seterr(divide='ignore')  # Disable numpy warning messages for dividing by 0

    read_list = [x.split('\t')[4:6] for x in open(inN).read().split('\n') if x!='' and x[0]!='#']
    read_list = [(int(x[1]), int(x[0])+int(x[1])) for x in read_list]
    
    # print >>sys.stderr, 'Enumerating log likelihoods'
    
    from numpy import log  # Faster function call
    from scipy.stats import binom
    from scipy.stats import poisson

    binom_pmf = binom.pmf
    poisson_pmf = poisson.pmf

    def compute_likelihood(alt_reads, num_reads, normal_freq, normal_alpha, cancer_freq, cancer_alpha, error_rate):
        alt_freq = round(normal_freq*normal_alpha + cancer_freq*cancer_alpha, 4)

        def pmf_mod(x, N, p):
            if x < 0 or x > N:
                return 0
            elif (x==0 and p==0) or (x==N and p==1):
                return 1
            else:
                return binom_pmf(x,N,p)

        mu = num_reads * error_rate  # Poisson parameter lambda

        x = 0
        x += sum([pmf_mod(alt_reads + k, num_reads, alt_freq)*poisson_pmf(k, mu) for k in range(0,min(2, num_reads-alt_reads)+1)])
        x += sum([pmf_mod(alt_reads - k, num_reads, alt_freq)*poisson_pmf(k, mu) for k in range(1,min(2, alt_reads)+1)])
        return x

    read_list = [(alt_reads, num_reads) for alt_reads, num_reads in read_list if alt_reads > cov_cutoff and (num_reads - alt_reads) > cov_cutoff and num_reads > cov_cutoff]
            
    f = open(inN.split('.txt')[0]+'.ll', 'w')

    ll_list = []
    for normal_alpha, cancer_alpha in parL:

        ll = sum([\
                log(sum([freq_prob*compute_likelihood(alt_reads, num_reads, normal_freq, normal_alpha, cancer_freq, cancer_alpha, error_rate)\
                             for normal_freq, cancer_freq, freq_prob in genL]))\
                    for alt_reads, num_reads in read_list if alt_reads!=0 and alt_reads!=num_reads])

        print >>f, str(ll)+'\t'+str(normal_alpha)+'\t'+str(cancer_alpha)
        f.flush()

        ll_list.append(ll)

    f.close()

    best_ll, best_par = max(zip(ll_list, parL), key=lambda x: x[0])                                

    ### Determine genotypes
    print >>outF, '#log-likelihood\t{0}\n#population frequencies\t{1}'.format(best_ll,'\t'.join([ str(x) for x in best_par ]))
    
    # Re-read the input tsv file (before imposing coverage cutoff)
    site_list = [x.split('\t') for x in open(inN).read().split('\n') if x!='' and x[0]!='#']

    parser.add_argument('infile', help='Input tsv file. Columns should be: chrom, position, ref base, alt base, number of reads supporting reference, number of reads supporting alternate.')
    
    for chrom, pos, t1, t2, ref_reads, alt_reads in site_list:
        alt_reads = int(alt_reads)
        num_reads = int(ref_reads) + alt_reads
        
        ll_list = [freq_prob*compute_likelihood(alt_reads, num_reads, normal_freq, best_par[0], cancer_freq, best_par[1], error_rate) for normal_freq, cancer_freq, freq_prob in genL]
        gen_idx, ll = max(enumerate(ll_list), key=lambda x:x[1])
        
        normal_gen, cancer_gen = genotypes[gen_idx]
        
        print >>outF, '\t'.join([chrom,pos,normal_gen,cancer_gen])

if __name__ == '__main__':
    main()
