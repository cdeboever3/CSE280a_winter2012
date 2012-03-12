import sys

print >>sys.stdout, """#!/bin/bash

# Compares coverage vs. alpha"""

covL = [25,50,75]
alphaL = [0.01,0.05,0.1,0.25]

for c in covL:
    for a in alphaL:

        tN = 'c_{0}_a_{1}.tsv'.format(c,a)
        rN = 'c_{0}_a_{1}_truth.txt'.format(c,a)
        gN = 'c_{0}_a_{1}_genotypes.txt'.format(c,a)

        print >>sys.stdout, '\n\n# coverage: {0} alpha: {1}\n\n'.format(c,a)

        print >>sys.stdout, 'python generate_sample_tsv_complex.tsv {0} -o {1} -ac {2} -sl 1000000 -sL {3} {4}\n'.format(rN,tN,c,a,1-a)
        
        print >>sys.stdout, 'python mixed_variant_calling.py {1} -er 0.001 > {0}\nwait\n'.format(gN,tN)
