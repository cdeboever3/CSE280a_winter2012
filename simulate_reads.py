#!/usr/bin/python

import argparse, subprocess, random, os, tempfile, threading, multiprocessing, shutil
import numpy as np
from sample_genome import weighted_random, format_fasta

def simulate_reads(n, ref_file, tmp_dir):
    tmp_dir = tempfile.gettempdir()  # Temp directory ('/tmp' on UNIX)
    tmp1 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write first reads
    tmp2 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write second reads

    # Sample reads with wgsim
    cmd = 'wgsim -N %s %s %s %s ' % (n, ref_file, tmp1, tmp2)
    p = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'w'),stderr=subprocess.STDOUT)
    p.wait()

    return tmp1, tmp2

def simulate_reads_star(x):
    return simulate_reads(*x)

def main():
    parser = argparse.ArgumentParser(description='Simulate a mixture of reads from multiple genomes using wgsim.')
    parser.add_argument('--genomes', nargs='+', required=True, help='FASTA files of genomes')
    parser.add_argument('--alpha', type=float, nargs='+', required=True, help='Mixing rates of reads (must sum up to 1). The i-th rate corresponds to the i-th genome listed for --genomes')
    parser.add_argument('--reads', required=True, type=float, help='The number of reads to sample')                      
    parser.add_argument('--e', type=float, default=0.02, help='Sequencing error rate when generating reads (Default: 0.02)')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--reads1', required=True, help='Output file for the first read of every paired-end.')
    parser.add_argument('--reads2', required=True, help='Output file for the second read of every paired-end.')
    parser.add_argument('--p', default=1, type=int, help='# of parallel processes')
    parser.add_argument('--tmp-dir', help='Temporary directory to write files (Default: create and use a randomly named directory within the current directory.)')
    args = parser.parse_args()

    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    try: args.reads = int(args.reads)
    except: raise Exception('Input an integer for --reads')
    assert args.p >=1, 'Enter positive number of threads'
    tmp_dir = tempfile.mkdtemp(dir='.') if (args.tmp_dir is None) else args.tmp_dir    
    if not os.path.isdir(tmp_dir): os.makedirs(tmp_dir)

    # Delete read files if they exist
    if os.path.isfile(args.reads1): os.remove(args.reads1)
    if os.path.isfile(args.reads2): os.remove(args.reads2)

    for f, a in zip(args.genomes, args.alpha):

        n = int(a * args.reads)  # Number of reads

        # Partition number of reads into roughly equal parts
        n_parts = [n / args.p for i in range(args.p)]
        n_parts[-1] += n % args.p
        
        # Run pool
        pool = multiprocessing.Pool(processes=args.p)
        tmp_files = pool.map(simulate_reads_star, [(x, f, tmp_dir) for x in n_parts])
        
        for (tmp1, tmp2), n_i in zip(tmp_files, n_parts):
            # Append reads to output files
            subprocess.call('cat %s | head -%s >> %s' % (tmp1, n_i*4, args.reads1), shell=True)
            subprocess.call('cat %s | head -%s >> %s' % (tmp2, n_i*4, args.reads2), shell=True)
        
    # Remove temp dir and temp files
    shutil.rmtree(tmp_dir)

if __name__=='__main__':
    main()   
