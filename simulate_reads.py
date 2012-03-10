#!/usr/bin/python

import argparse, subprocess, random, os, tempfile, threading, multiprocessing, shutil, sys, pdb
import numpy as np
from sample_genome import weighted_random, format_fasta

def simulate_reads(n, error, ref_file, tmp_dir, proc_index, debug):
    # tmp_dir = tempfile.gettempdir()  # Temp directory ('/tmp' on UNIX)
    tmp1 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write first reads
    tmp2 = tempfile.mkstemp('.fq', dir=tmp_dir)[1]  # File to write second reads

    codeD = os.path.dirname(os.path.realpath(__file__)) # directory containing this script and other needed scripts

    # Sample reads with wgsim
    cmd = '%s -N %s -e %s %s %s %s ' % (os.path.join(codeD, 'wgsim'), n, error, ref_file, tmp1, tmp2)
    if debug: print >> sys.stderr, "Process %i: Starting wgsim: %s\n" %(proc_index, cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=subprocess.PIPE)    
    if debug: print >> sys.stderr, 'Process %i: Waiting for wgsim to finish...' % proc_index
    p.wait()
    if debug: print >> sys.stderr, 'Process %i: wgsim finished.' % proc_index
    if debug: print >> sys.stderr, 'Process %i: wgsim stderr:\n' % proc_index, p.communicate()[1]    
    if p.returncode != 0:
        print >> sys.stderr, "Process %i: wgsim did not run properly (it's shell exit code was not 0)." % proc_index
        assert p.returncode==0

    return tmp1, tmp2

def simulate_reads_star(x):
    return simulate_reads(*x)

def main():
    parser = argparse.ArgumentParser(description='Simulate a mixture of reads from multiple genomes using wgsim.')
    parser.add_argument('--genomes', nargs='+', required=True, help='FASTA files of genomes')
    parser.add_argument('--alpha', type=float, nargs='+', required=True, help='Mixing rates of reads (must sum up to 1). The i-th rate corresponds to the i-th genome listed for --genomes')
    parser.add_argument('--reads', required=True, type=float, help='Number of reads to sample')
    parser.add_argument('--e', type=float, default=0.02, help='Sequencing error rate when generating reads (default: 0.02)')
    parser.add_argument('--reads1', required=True, help='Output file for the first read of every paired-end.')
    parser.add_argument('--reads2', required=True, help='Output file for the second read of every paired-end.')
    parser.add_argument('--p', default=1, type=int, help='# of parallel processes (default: 1)')
    parser.add_argument('--tmp-dir', default='.', help='Temporary directory to write intermediate read files (default: create and use a randomly named directory within the current directory.)')
    parser.add_argument('--ram-disk', help='A RAM-mounted directory to write temporary copies of the genomes and final reads.')
    parser.add_argument('--debug', action='store_true', help='Print debug statements to stdout')
    args = parser.parse_args()

    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    try: args.reads = int(args.reads)
    except: raise Exception('Input an integer for --reads')
    assert args.p >=1, 'Enter positive number of threads'

    # Make temp dir for intermediate read files
    if args.debug:
        print >>sys.stderr, 'tmp_dir: {0}'.format(args.tmp_dir)
    tmp_dir = tempfile.mkdtemp(dir=('.' if (args.tmp_dir is None) else args.tmp_dir))
    if not os.path.isdir(tmp_dir): os.makedirs(tmp_dir)

    if args.ram_disk is None:
        genomes = args.genomes
        # reads1 = args.reads1
        # reads2 = args.reads2
    else:
        # Copy genomes to RAM DISK
        if args.debug: print >> sys.stderr, "Copying genomes over to RAM disk.."
        ram_genomes = [os.path.join(args.ram_disk, os.path.basename(x)) for x in args.genomes]
        for g, ram_g in zip(args.genomes, ram_genomes):  shutil.copy(g, ram_g)
        genomes = ram_genomes
        
        # Temporary files to consolidate intermediate read files
        # reads1 = os.path.join(args.ram_disk, os.path.basename(args.reads1))
        # reads2 = os.path.join(args.ram_disk, os.path.basename(args.reads2))
        
    reads1 = args.reads1
    reads2 = args.reads2

    # Delete read files if they exist
    if os.path.isfile(args.reads1): os.remove(args.reads1)
    if os.path.isfile(args.reads2): os.remove(args.reads2)

    for g, a in zip(genomes, args.alpha):

        n = int(a * args.reads)  # Number of reads

        # Partition number of reads into roughly equal parts
        n_parts = [n / args.p for i in range(args.p)]
        n_parts[-1] += n % args.p
        
        # Run process pool
        pool = multiprocessing.Pool(processes=args.p)
        tmp_files = pool.map(simulate_reads_star, [(x, args.e, g, tmp_dir, i, args.debug) for i, x in enumerate(n_parts)])

        if args.debug: 
            sys.stderr.flush()
            print >> sys.stderr, "All processes finished running wgsim.\nConcatenating intermediate read files.."

        # Concatenate intermediate read files        
        for (tmp1, tmp2), n_i in zip(tmp_files, n_parts):
            # Lines are filtered with head because sometimes wgsim
            # produces 1 extra paired end for unknown reasons
            p = subprocess.Popen('cat %s | head -%s >> %s' % (tmp1, n_i*4, reads1), shell=True)
            p.wait()
            assert p.returncode==0, 'Concatenating intermediate read file %s to final read file %s failed' % (tmp1, reads1)
            p = subprocess.Popen('cat %s | head -%s >> %s' % (tmp2, n_i*4, reads2), shell=True)
            p.wait()            
            assert p.returncode==0, 'Concatenating intermediate read file %s to final read file %s failed' % (tmp2, reads2)
    
    if args.debug: print >> sys.stderr, "Removing intermediate read files.."
    shutil.rmtree(tmp_dir)  # Remove temp dir and temp files
    
    if args.ram_disk is not None:
        if args.debug: print >> sys.stderr, "Removing copies of genomes in RAM folder.."
        # Remove RAM copies of genomes and read files
        for g in genomes: os.remove(g)
        # if args.debug: print >> sys.stderr, "Moving final read file from RAM to desired folder.."
        # shutil.move(reads1, args.reads1)
        # shutil.move(reads2, args.reads2)


if __name__=='__main__':
    main()   
