#!/usr/bin/python

import argparse, subprocess, random, os, tempfile
import numpy as np
from sample_genome import weighted_random, format_fasta

def main():
    parser = argparse.ArgumentParser(description='Simulate a mixture of reads from multiple genomes.')
    parser.add_argument('--genomes', nargs='+', required=True, help='FASTA files of genomes')
    parser.add_argument('--alpha', type=float, nargs='+', required=True, help='Mixing rates of reads (must sum up to 1). The i-th rate corresponds to the i-th genome listed for --genomes')
    parser.add_argument('--reads', required=True, type=float, help='The number of reads to sample')                      
    parser.add_argument('--e', type=float, default=0.02, help='Sequencing error rate when generating reads (Default: 0.02)')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--reads1', required=True, help='Output file for the first read of every paired-end.')
    parser.add_argument('--reads2', required=True, help='Output file for the second read of every paired-end.')
    args = parser.parse_args()

    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    try: args.reads = int(args.reads)
    except: raise Exception('Input an integer for --reads')
        
    if os.path.isfile(args.reads1): os.remove(args.reads1) 
    if os.path.isfile(args.reads2): os.remove(args.reads2) 

    for f, a in zip(args.genomes, args.alpha):

        n = int(a * args.reads)  # Number of reads
        
        tmp_dir = tempfile.gettempdir()  # Temp directory ('/tmp' on UNIX)
        tmp1 = os.path.join(tmp_dir, os.path.basename(f + '.1.fq'))  # File to write first reads
        tmp2 = os.path.join(tmp_dir, os.path.basename(f + '.2.fq'))  # File to write second reads

        # Sample reads with wgsim
        cmd = 'wgsim -N %s %s %s %s ' % (n, f, tmp1, tmp2)
        p = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'w'),stderr=subprocess.STDOUT)
        p.wait()

        # Append reads to output file
        subprocess.call('cat %s | head -%s >> %s' % (tmp1, n*4, args.reads1), shell=True)
        subprocess.call('cat %s | head -%s >> %s' % (tmp2, n*4, args.reads2), shell=True)

if __name__=='__main__':
    main()   
