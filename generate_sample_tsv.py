# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, random
import numpy as np
from numpy.random import binomial
from random import gauss

### helper functions ###

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

if __name__ == '__main__':

    ploidyL = [2,2]
    subL = [0.25,0.75]
    num_sites = 200
    num_noise_sites = 10
    error_rate = 0.01
    avg_cov = 50
    std_cov = 10

    parser = argparse.ArgumentParser(description='This script makes a sample input file for the variant calling script by simulating read counts using a binomial distribution and coverage with a Gaussian.')
    parser.add_argument('-o', nargs='?', type=argparse.FileType('w'),default=sys.stdout, help='Output file. Default: standard out')
    parser.add_argument('-pL', default=ploidyL, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-sL', default=subL, nargs='+', help='List of subpopulation frequencies. Default: {0}'.format(' '.join([ str(x) for x in subL ])))
    parser.add_argument('-ns', default=num_sites, type=int, help='Number of variable sites. Default: {0}'.format(num_sites))
    parser.add_argument('-nn', default=num_noise_sites, type=int, help='Number of variable sites that are noise. Default: {0}'.format(num_noise_sites))
    parser.add_argument('-er', default=error_rate, type=int, help='Error rate for noise sites. Default: {0}'.format(error_rate))
    parser.add_argument('-ac', default=avg_cov, type=int, help='Average coverage. Default: {0}'.format(avg_cov))
    parser.add_argument('-sc', default=std_cov, type=int, help='Standard deviation of coverage. Default: {0}'.format(std_cov))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
   
    outF            = args.o
    ploidyL         = args.pL
    subL            = args.sL
    num_sites       = args.ns
    num_noise_sites = args.nn
    error_rate      = args.er
    avg_cov         = args.ac
    std_cov         = args.sc
    debug           = args.d

    if len(ploidyL) > 2:
        print >>sys.stderr, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    exp_freqL = generate_possible_freqL(ploidyL,subL)[1:] # we don't want 0 frequency

    for iii in range(num_sites-num_noise_sites):
        cov = int(round(gauss(avg_cov,std_cov)))
        alt_freq = random.choice(exp_freqL)
        alt_reads = binomial(cov,alt_freq)

        print >>outF, '{0}\t{1}'.format(cov-alt_reads,alt_reads)
    for iii in range(num_noise_sites):
        cov = int(round(gauss(avg_cov,std_cov)))
        alt_reads = binomial(cov,error_rate)

        print >>outF, '{0}\t{1}'.format(cov-alt_reads,alt_reads)
