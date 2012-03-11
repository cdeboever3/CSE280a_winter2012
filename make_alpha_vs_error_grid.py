import sys

print >>sys.stdout, """#!/bin/bash

# Compares alpha vs. error

#PBS -N alpha_vs_error_grid
#PBS -l nodes=1:ppn=8
#PBS -o /raid/development/cdeboever/CSE280a/alpha_vs_error_grid/alpha_vs_error_grid.out
#PBS -e /raid/development/cdeboever/CSE280a/alpha_vs_error_grid/alpha_vs_error_grid.err

mkdir -p /scratch/alpha_vs_error_grid

cd /scratch/alpha_vs_error_grid

rsync -avz /raid/development/cdeboever/CSE280a/c1.fa .
rsync -avz /raid/development/cdeboever/CSE280a/n1.fa .

wait

"""

alphaL = [0.01]
errorL = [0.001]
# alphaL = [0.01,0.05,0.1,0.25]
# errorL = [0.001, 0.01, 0.02]

for a in alphaL:
    for e in errorL:

        tN = 'a_{0}_e_{1}.tsv'.format(a,e)
        gN = 'a_{0}_e_{1}_genotypes.txt'.format(a,e)

        print >>sys.stdout, '\n\n# alpha: {0} error: {1}\n\n'.format(a,e)

        print >>sys.stdout, 'python /raid/development/cdeboever/CSE280a/CSE280a_winter2012/simulate_reads.py --genomes n1.fa c1.fa --alpha {0} {1} --reads 11428571 --e {2} --reads1 a01r11428571e001.R1.fastq --reads2 a01r11428571e001.R2.fastq --ram-disk /dev/shm --tmp-dir /scratch/alpha_vs_error_grid\nwait\n'.format(a,1-a,e)
        
        print >>sys.stdout, 'python /raid/development/cdeboever/CSE280a/CSE280a_winter2012/make_variant_calling_input.py a01r11428571e001.R1.fastq a01r11428571e001.R2.fastq -o {0}\nwait\n'.format(tN)
        
        print >>sys.stdout, 'python /raid/development/cdeboever/CSE280a/CSE280a_winter2012/mixed_variant_calling.py {2} -pL 2 2 -er {0} > {1}\nwait\n'.format(e,gN,tN)

print >>sys.stdout,"""
rsync -avz *tsv /raid/development/cdeboever/CSE280a/alpha_vs_error_grid
rsync -avz *genotypes.txt /raid/development/cdeboever/CSE280a/alpha_vs_error_grid

wait
 
rm -r /scratch/alpha_vs_error_grid
"""
