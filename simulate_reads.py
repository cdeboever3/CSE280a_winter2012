#!/usr/bin/python

import argparse, subprocess, random, os, tempfile, threading, multiprocessing, shutil
import numpy as np
from sample_genome import weighted_random, format_fasta

def simulate_reads(n, ref_file, tmp_dir, ram_disk=None):
    tmp_dir = tempfile.gettempdir()  # Temp directory ('/tmp' on UNIX)
    tmp1 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write first reads
    tmp2 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write second reads

    codeD = os.path.dirname(os.path.realpath(__file__)) # directory containing this script and other needed scripts

    # Sample reads with wgsim
    cmd = '%s -N %s %s %s %s ' % (os.path.join(codeD, 'wgsim'), n, ref_file, tmp1, tmp2)
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
    parser.add_argument('--tmp-dir', default='.', help='Temporary directory to write files (Default: create and use a randomly named directory within the current directory.)')
    parser.add_argument('--ram-disk', help='Temporarily copy the genome files and write read files to this RAM-mounted directory for faster I/O.  If not specified, then genome files are read from their original location.')
    args = parser.parse_args()

    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    try: args.reads = int(args.reads)
    except: raise Exception('Input an integer for --reads')
    assert args.p >=1, 'Enter positive number of threads'
    tmp_dir = tempfile.mkdtemp(dir=('.' if (args.tmp_dir is None) else args.tmp_dir))
    if not os.path.isdir(tmp_dir): os.makedirs(tmp_dir)

    print tmp_dir

    # Copy genomes to RAM DISK
    if args.ram_disk is not None:
        ram_genomes = [os.path.join(args.ram_disk, os.path.basename(x)) for x in args.genomes]
        for g, ram_g in zip(args.genomes, ram_genomes):
            shutil.copy(g, ram_g)
        genomes = ram_genomes
        reads1 = os.path.join(args.ram_disk, os.path.basename(args.reads1))
        reads2 = os.path.join(args.ram_disk, os.path.basename(args.reads2))
    else:
        genomes = args.genomes
        reads1 = args.reads1
        reads2 = args.reads2

    # Delete read files if they exist
    if os.path.isfile(args.reads1): os.remove(args.reads1)
    if os.path.isfile(args.reads2): os.remove(args.reads2)

    for f, a in zip(genomes, args.alpha):

        n = int(a * args.reads)  # Number of reads

        # Partition number of reads into roughly equal parts
        n_parts = [n / args.p for i in range(args.p)]
        n_parts[-1] += n % args.p
        
        # Run pool
        pool = multiprocessing.Pool(processes=args.p)
        tmp_files = pool.map(simulate_reads_star, [(x, f, tmp_dir, args.ram_disk) for x in n_parts])
        
        for (tmp1, tmp2), n_i in zip(tmp_files, n_parts):
            # Append reads to output files
            subprocess.call('cat %s | head -%s >> %s' % (tmp1, n_i*4, reads1), shell=True)
            subprocess.call('cat %s | head -%s >> %s' % (tmp2, n_i*4, reads2), shell=True)
        
    # Remove temp dir and temp files
    shutil.rmtree(tmp_dir)
    # Remove RAM copies of genomes and read files
    if args.ram_disk is not None:
        for g in genomes: os.remove(g)
        shutil.move(reads1, args.reads1)
        shutil.move(reads2, args.reads2)

if __name__=='__main__':
    main()   
