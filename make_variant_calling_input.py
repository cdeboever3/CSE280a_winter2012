# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, glob, os
import numpy as np

### helper functions ###

if __name__ == '__main__':

    ### magic variables ###
    # these variables can be set at the command line as well
    keep_bam = False
    keep_pileup = False
    reference = '/raid/references-and-indexes/hg19/bwa_new/hg19.fa'
    seedlen = 35

    ### gather command line arguments ###
    parser = argparse.ArgumentParser(description='This script takes two fastq files representing R1 and R2 reads, aligns them to hg19, creates a pileup file, and parses the pileup file to make the tsv file needed for mixed_variant_calling.py.')
    parser.add_argument('R1_reads', help='Input fastq file with R1 reads.')
    parser.add_argument('R2_reads', default=None, help='Input fastq file with R2 reads.')
    parser.add_argument('-kb', action='store_true', default=keep_bam, help='Keep bam file. Default: {0}'.format(keep_bam))
    parser.add_argument('-kp', action='store_true', default=keep_pileup, help='Keep pileup file. Default: {0}'.format(keep_pileup))
    parser.add_argument('-r', metavar='reference', help='Reference for BWA. Default {0}'.format(reference))
    parser.add_argument('-s', metavar='seedlen', type=int, help='Seed length for BWA. Default: {0}'.format(seedlen))
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
    
    R1N         = args.R1_reads
    R2N         = args.R2_reads
    keep_bam    = args.kb
    keep_pileup = args.kp
    reference   = args.r
    seedlen     = args.s
    debug       = args.d

    R1N = os.path.realpath(R1N) # get the input file path
    R2N = os.path.realpath(R2N) # get the input file path

    ### make file names ###
    prefix = '.'.join(R1N.split('.')[1:])
    bamN = prefix + '.bam'
    pileupN = prefix + '.pileup'
    outputN = prefix + '_counts.tsv'

    ### align with BWA ###

    p = subprocess.Popen('bwa sampe {0} <(bwa aln -t 4 -l {1} {0} {2}) <(bwa aln -t 4 -l {1} {0} {3}) {2} {3} | samtools view -Sb - | samtools sort -o temp_read_sort > {4}'.format(reference, seedlen, R1N, R2N, outputN),shell=True)
    sts = os.waitpid(p.pid, 0)[1] # wait for process to finish

    ### make pileup file ###

    # p = subprocess.Popen('samtools view -c {0}{1}.bam'.format(basename,),shell=True)
    # sts = os.waitpid(p.pid, 0)[1] # wait for process to finish

    # ### make input file for mixed_variant_calling.py ###

    # p = subprocess.Popen('samtools view -c {0}{1}.bam'.format(basename,),shell=True)
    # sts = os.waitpid(p.pid, 0)[1] # wait for process to finish

    # ### remove bam and pileup files ###

    # if not keep_bam:
    #     os.remove(bamN)
    # if not keep_pileup:
    #     os.remove(pileupN)
