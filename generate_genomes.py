#!/usr/bin/python

import argparse, subprocess, random, os
import numpy as np
from sample_genome import weighted_random, format_fasta

def generate_allele_freq(alphabet, distribution):
    """Generate allele frequences.
    @alphabet : List of alleles
    @distribution : Distribution to sample frequencies."""

    if distribution.upper()=='UNIFORM':
        freq_list = [random.random() for i in range(len(alphabet))]
    else:
        raise Exception("Distribution %s isn't supported" % distribution)
    
    # Normalize to sum to 1
    tot_freq = sum(freq_list)
    freq_list = [float(x)/tot_freq for x in freq_list]
    
    return freq_list

def mutate_seq(seq, mu, pl):
    """
    Replicate a sequence <seq> into <pl> copies.  Mutate the same
    position in all copies at rate mu.
    
    @seq : sequence to mutate
    @mu  : mutation rate
    @pl  : ploidy number
    """
    n = len(seq)         # Sequence length
    alphabet = list(set(seq))  # Sequence alphabet: A,T,C,G in DNA
    
    # Number of mutations, drawn from binomial dist
    num_mut = np.random.binomial(n, mu)

    # Randomly select mutation positions uniformly at random
    pos_list = sorted(random.sample(range(n), num_mut))
    
    seq_list = [list(seq) for i in range(pl)]  # Replicate sequence

    # Mutate positions:
    for pos in pos_list:
        
        # Generate new alleles for position
        orig_alleles = [s[pos] for s in seq_list]
        freq_list = generate_allele_freq(alphabet, 'UNIFORM')
        while True:
            new_alleles = [weighted_random(alphabet, freq_list) for i in range(pl)]
                        
            # Repeat if every new allele matches the original allele
            if not all([x==y for x,y in zip(new_alleles, orig_alleles)]): break
        
        if debug: print pos+1, orig_alleles, new_alleles, freq_list

        # Mutate
        for a, seq in zip(new_alleles, seq_list): seq[pos] = a

    seq_list = [''.join(s) for s in seq_list]
    return seq_list

def main():
    parser = argparse.ArgumentParser(description='Generate normal or cancer genomes by mutating a reference sequence')
    parser.add_argument('--ref', required=True, help='FASTA file of a contiguous reference sequence')
    parser.add_argument('--k', default=['n1, c1'], nargs='+', required=True, help="List of genome types to generate.  Each genome type is denoted by a two-letter word where the first letter is either 'n' (normal) or 'c' (cancer).  The second letter is a positive integer indicating which individual the genome comes from.  For example, 'n1 c1 n2' means generating three genomes where the first two genomes are a normal and cancer genome from the same individual, and where the third genome is a normal one taken from a different individual.  As another example, 'n1 n2 n3 n4' means generating four normal genomes, each of which comes from a different individual.")
    parser.add_argument('--pl', required=True, type=int, help='Genome ploidy')
    parser.add_argument('--germ-mu', type=float, required=True, help='Germline mutation rate, i.e. the probability that any haploid in a normal genome differs from the reference at a particular site.')
    parser.add_argument('--som-mu', type=float, required=True, help='Somatic x mutation rate, i.e. the probability that any haploid in a cancerous genome differs from the reference at a particular site.')
    parser.add_argument('--alpha', type=float, nargs='+', required=True, help='Mixing rates of reads.  ')
    parser.add_argument('--reads', required=True, type=float, help='The number of reads to sample')                      
    parser.add_argument('--e', type=float, default=0.02, help='Sequencing error rate when generating reads (Default: 0.02)')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--reads1', required=True, help='Output file for the first read of every paired-end.')
    parser.add_argument('--reads2', required=True, help='Output file for the second read of every paired-end.')
    args = parser.parse_args()

    global debug
    debug = args.debug

    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    assert len(args.alpha)==args.pl, "Enter same number mixing rates as ploidy"
    for x in args.k:
        assert x[0] in ['n', 'c'], "Bad format for --k.  Genome type must either be 'n' for normal or 'c' for cancer"        
        try: int(x[1])
        except ValueError: raise Exception("Bad format for --k. Use positive integer to denote individual")
        assert len(x)==2, "Bad format for --k.  Enter two-letter genome type descriptors."
    try: args.reads = int(args.reads)
    except: raise Exception('Input an integer for --reads')

    # Read reference genomic region
    header = open(args.ref).read().split('\n')[0]
    assert '>' in header, 'First line of FASTA file needs to be the header'
    header = header[1:]  # Remove '>'
    ref_seq = ''.join(open(args.ref).read().split('\n')[1:])  
    individuals = set([x[1] for x in args.k])

    # Dict: i --> normal genome of individual i
    normal_dict = dict()

    # Sort genome types by individual indices (used later for deleting
    # dictionary entries)
    args.k.sort(key=lambda x: x[1])

    # Construct genomes
    i = None
    for d in args.k:
        
        # Delete dictionary entry to save memory 
        if (i is not None) and d[1]!=i: del normal_dict[i]

        t = d[0]  # Type of genome: normal or cancer
        i = d[1]  # Index of individual        
        
        if debug: print ('Cancer' if t=='c' else 'Normal') + ' genome for individual %s' %i
        # Construct normal genome of individual if it's requested or
        # needed to construct cancer genome
        if t=='n' or (t=='c' and not normal_dict.has_key(i)):  
            seq_list = mutate_seq(ref_seq, args.germ_mu, args.pl)
            normal_dict[i] = seq_list

        # Construct cancer genome by mutating normal genome of the same individual        
        if t=='c':
            seq_list = [mutate_seq(s, args.som_mu, 1)[0] for s in normal_dict[i]]

        if debug:
            for j in range(args.pl): print 'Mutations in Haploid %s: %s' % (j, sum([x!=y for x,y in zip(ref_seq, seq_list[j])]))

        # Write genome to FASTA file
        fasta = '\n'.join([format_fasta('%s | haploid %s' % (header,j+1), s) for j, s in enumerate(seq_list)])
        open(d+'.fa', 'w').write(fasta)

        
    if os.path.isfile(args.reads1): os.remove(args.reads1) 
    if os.path.isfile(args.reads2): os.remove(args.reads2) 

    for d, a in zip(args.k, args.alpha):

        n = int(a * args.reads)  # Number of reads
        
        f = d + '.fa'  # Genome file
        tmp1 = d + '.1.fq'  # File to write first reads
        tmp2 = d + '.2.fq'  # File to write second reads

        # Sample reads with wgsim
        cmd = 'wgsim -N %s %s %s %s ' % (n, f, tmp1, tmp2)
        p = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'w'),stderr=subprocess.STDOUT)
        p.wait()

        # Append first reads to output file
        subprocess.call('cat %s | head -%s >> %s' % (tmp1, n*4, args.reads1), shell=True)
        subprocess.call('cat %s | head -%s >> %s' % (tmp2, n*4, args.reads2), shell=True)

if __name__=='__main__':
    main()   
