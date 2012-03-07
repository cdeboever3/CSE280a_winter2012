#!/usr/bin/python
import argparse, subprocess, random
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
    freq_list = freq_list / len(freq_list)  # Normalize to sum to 1
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
    alphabet = set(seq)  # Sequence alphabet: A,T,C,G in DNA
    
    # Number of mutations, drawn from binomial dist
    num_mut = np.random.binomial(n, args.mu)
    # Randomly select mutation positions uniformly at random
    pos_list = random.sample(range(n), num_mut)
    
    seq_list = [seq for seq in range(pl)]  # Replicate sequence
    
    # Mutate positions:
    for pos in pos_list:
        
        # Generate new alleles for position
        while True:
            freq_list = generate_allele_freq(alphabet, 'UNIFORM')
            new_alleles = [weighted_random(alphabet, freq_prob) for i in range(pl)]
            
            # Repeat if every new allele matches the original allele
            orig_alleles = [s[pos] for s in seq_list]
            if not all([x==y for x,y in zip(new_alleles, orig_alleles)]): break
        
        # Mutate
        for a, seq in zip(new_alleles, seq_list): seq[pos] = a
    
    return seq_list
def main():
    parser = argparse.ArgumentParser(description='Generate normal or cancer genomes by mutating a reference sequence')
    parser.add_argument('--ref', required=True, help='FASTA file of a contiguous reference sequence')
    parser.add_argument('--k', default=['n1, c1'], help="List of genome types to generate.  Each genome type is denoted by a two-letter word where the first letter is either 'n' (normal) or 'c' (cancer).  The second letter is a positive integer indicating which individual the genome comes from.  For example 'n1 c1 n2' means that three genomes are generated where the first two genomes are a normal and cancer genome, respectively, from the same individual, and where the third genome is a normal one taken from a different individual.  As another example 'n1 n2 n3 n4' means generating four normal genomes, each of which comes from a different individual.")
    parser.add_argument('--pl', required=True, type=int, help='Ploidy number')
    parser.add_argument('--germ-mu', type=float, required=True, help='Germline mutation rate, i.e. the probability that any haploid in a normal genome differs from the reference at a particular site.')
    parser.add_argument('--som-mu', type=float, required=True, help='Somatic mutation rate, i.e. the probability that any haploid in a cancerous genome differs from the reference at a particular site.')
    parser.add_argument('--alpha', type=float, nargs='+', required=True, help='Mixing rates of reads.  ')
    parser.add_argument('--reads', type=int, required=True, help='The number of reads to sample')                      
    parser.add_argument('--e', type=float, default=0.02, help='Sequencing error rate when generating reads (Default: 0.02)')
    parser.add_argument('--output', required=True, help='Output file')
    args = parser.parse_args()
    
    # Check params
    assert sum(args.alpha)==1, "Mixing rates don't add up to 1"
    assert len(args.alpha)==args.pl, "Enter same number mixing rates as ploidy"
    assert set([x[0] for x in args.k])==set(['n', 'c']), "Genome types must either be 'n' for normal or 'c' for cancer"
    for x in args.k:
        assert x[0] in ['n', 'c'], "Bad format for --k.  Genome type must either be 'n' for normal or 'c' for cancer"
        
        try: int(x[1])
        except ValueError: raise Exception("Bad format for --k. Use positive integer to denote individual")
        
        assert len(x)==2, "Bad format for --k.  Enter two-letter genome type descriptors."
    # Read reference genomic region
    header = open(args.ref).read().split('\n')[0]
    assert '>' in header, 'First line of FASTA file needs to be the header'
    seq = ''.join(open(args.ref).read().split('\n')[1:])  
    individuals = set([x[1] for x in args.k])
    # Dict: i --> normal genome of individual i
    normal_dict = dict()
    for i in individuals:
        # Construct normal genome of individual with ploidy pl
        seq_list = mutate_seq(seq, args.germ_mu, args.pl)
        # Write to FASTA
        individual_file = 'n%i.fa' % i
        fasta = '\n'.join([format_fasta(header, seq_list[i]) for i in range(seq)])
        open(individual_file, 'w').write(fasta)
        
        normal_dict[i] = seq_list
        
    for i in [x[1] for x in args.k if x[0]=='c']:
        
        # Construct a cancer genome from the normal genom of the same individual
        cancer_file = 'c%i.fa' % i
        seq_list = [mutate_seq(seq, args.som_mu, 1)[0] for seq in normal_dict[i]]
        fasta = '\n'.join([format_fasta(header, seq_list[j]) for j in range(seq)])
        open(cancer_file, 'w').write(fasta)
    for t, a in zip(args.k, args.alpha):
        n = a * args.reads         # Number of reads
        f = '%s%i' % (t[0], t[1])  # Genome file
        tmp1 = 'tmp1.fq'  # File to write first reads
        tmp2 = 'tmp2.fq'  # File to write second reads
        # Sample reads with wgsim
        p = subprocess.Popen('wgsim -N %s %s %s %s ' % (n, f, tmp1, tmp2), shell=True, stdout=open(os.devnull, 'w'))
        p.wait()
        # Extract the first reads
        subprocess.call('cat %s >> %s' % (tmp1, args.output))
if __name__=='__main__':
    main()   
