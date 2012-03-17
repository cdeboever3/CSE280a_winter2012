import sys

print >>sys.stdout, """#!/bin/bash

# Compares alpha vs. error"""

alphaL = [0.01,0.05,0.1,0.25]
errorL = [0.001, 0.01, 0.02]

for a in alphaL:
    for e in errorL:

        tN = 'a_{0}_e_{1}.tsv'.format(a,e)
        rN = 'a_{0}_e_{1}_truth.txt'.format(a,e)
        gN = 'a_{0}_e_{1}_genotypes.txt'.format(a,e)

        print >>sys.stdout, '\n\n# alpha: {0} error: {1}\n\n'.format(a,e)

        print >>sys.stdout, 'python generate_sample_tsv_complex.tsv {0} -o {1} -sl 1000000 -er {2} -sL {3} {4}\n'.format(rN,tN,e,a,1-a)
        
        print >>sys.stdout, 'python mixed_variant_calling.py {2} -er {0} > {1}\nwait\n'.format(e,gN,tN)
