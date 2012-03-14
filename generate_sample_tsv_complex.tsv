# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, random
import sys, pdb, random
import numpy as np
from numpy.random import binomial
from scipy.stats import poisson
# from mixed_variant_calling import generate_possible_freqL

### helper functions ###

if __name__ == '__main__':

    ploidyL = [2,2]
    subL = [0.25,0.75]
    seq_len = 100000000	
    germ_mu = 0.001
    soma_mu = 0.00001
    error_rate = 0.001
    avg_cov = 50

    parser = argparse.ArgumentParser(description='This script makes a sample input file for the variant calling script by simulating read counts using a binomial distribution and coverage with a Poisson.')
    parser.add_argument('to', help='Output file for truth.')
    parser.add_argument('-o', nargs='?', type=argparse.FileType('w'),default=sys.stdout, help='Output file. Default: standard out')
    parser.add_argument('-pL', default=ploidyL, type=int, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
    parser.add_argument('-sL', default=subL, type=float, nargs='+', help='List of subpopulation frequencies. Default: {0}'.format(' '.join([ str(x) for x in subL ])))
    parser.add_argument('-sl', default=seq_len, type=int, help='Length of sequence to simulate. Default: {0}'.format(seq_len))
    parser.add_argument('-gm', default=germ_mu, type=int, help='Germline mutation rate. Default: {0}'.format(germ_mu))
    parser.add_argument('-sm', default=soma_mu, type=int, help='Somatic mutation rate. Default: {0}'.format(soma_mu))
    parser.add_argument('-er', default=error_rate, type=float, help='Error rate for noise sites. Default: {0}'.format(error_rate))
    parser.add_argument('-ac', default=avg_cov, type=int, help='Average coverage. Default: {0}'.format(avg_cov))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()

    truthN                    = args.to
    outF            = args.o
    ploidyL         = args.pL
    subL            = args.sL
    seq_len			= args.sl
    germ_mu			= args.gm
    soma_mu			= args.sm
    error_rate      = args.er
    avg_cov         = args.ac
    debug           = args.d

    if len(ploidyL) > 2:
        print >>sys.stderr, 'Sorry, only two subpopulations are currently supported.'
        sys.exit(1)

    truthF = open(truthN,'w')
    
    print >>outF, '#chromosome\tposition (0-based)\tref base\talt base\tref coverage\talt coverage'
    print >>truthF, '#log-likelihood\t{0}\n#population frequencies\t{1}'.format(0,'\t'.join([ str(x) for x in subL ]))

    for pos in xrange(seq_len):
	alt_freq = 0
        germline = 'AA'
        cancer = 'AA'
        cov = int(round(poisson.rvs(avg_cov)))

        if cov==0: continue  # Continue if base isn't covered

	if random.random() < germ_mu: # germline mutation, hetero or homozygous
	    alt_freq += random.choice([0.5,1])
            if alt_freq == 1:
                germline = 'TT'
                cancer = 'TT'
            else:
                germline = 'AT'
                cancer = 'AT'
	elif random.random() < soma_mu: # if not germline, somatic (infinite sites)
	    alt_freq += subL[1]*0.5
            cancer = 'AT'
        alt_reads = binomial(cov,alt_freq)
	ref_reads = cov - alt_reads
	for i in xrange(cov):
#	    if random.random < error_rate:
	    if random.random() < error_rate:
		d = random.choice([-1,1])
		alt_reads += d
#		ref_reads += d
		ref_reads -= d
            # if alt_reads < 0:
	    #     alt_reads = abs(alt_reads)
	    #     ref_reads -= alt_reads
	    # if ref_reads < 0:
	    #     ref_reads = abs(ref_reads)
	    #     alt_reads -= ref_reads
        if alt_reads < 0:
            ref_reads += alt_reads
            alt_reads = 0
        if ref_reads < 0:
            alt_reads += ref_reads
            ref_reads = 0

#        if alt_reads != 0:
        if alt_reads != 0 and ref_reads != 0:
            print >>outF, 'chr1\t{2}\tA\tT\t{0}\t{1}'.format(ref_reads,alt_reads,pos)
            print >>truthF, '\t'.join(['chr1',str(pos),germline,cancer])
