#!/usr/bin/python

import argparse, subprocess, random, os
import numpy as np
from sample_genome import weighted_random, format_fasta

def generate_allele_freq(allele_space, distribution):
    """Generate allele frequences.
    @allele_space : List of alleles
    @distribution : Distribution to sample frequencies."""

    if distribution.upper()=='UNIFORM':
        freq_list = [random.random() for i in range(len(allele_space))]
    else:
        raise Exception("Distribution %s isn't supported" % distribution)
    
    # Normalize to sum to 1
    tot_freq = sum(freq_list)
    freq_list = [float(x)/tot_freq for x in freq_list]
    
    return freq_list

def mutate_seq(ref, mu, pl, distribution='UNIFORM'):
    """
    Generate a <pl> ploidy genome by replicating and mutating a
    reference sequence.  For each mutation site, allow only one
    possible alternate allele.  For example if the reference allele at
    a site is 'A', then the corresponding site on both haploids of a
    heterozygote diploid genome cannot be ('C','G') where there are
    two alternate alleles.  Allele frequencies are randomly chosen
    according to a specified distribution
    
    @ref : reference sequence
    @mu  : mutation rate, i.e. probability that any haploid at a
           particular site differs from the reference
    @pl  : genome ploidy

    Returns (seq_list, mutation_list) where
    @ seq_list : List of genome's haploids
    @ mutation_dict[i] : genotype of mutation at position i (0-based)
    """

    n = len(ref)               # Sequence length
    alphabet = list(set(ref))  # Sequence alphabet: A,T,C,G in DNA
    
    # Number of mutations (drawn from binomial dist)
    num_mut = np.random.binomial(n, mu)

    # List of positions to mutate (chosen uniformly at random)
    pos_list = sorted(random.sample(range(n), num_mut))
    
    # Replicate reference for ploidy genome
    seq_list = [list(ref) for i in range(pl)]
    
    # dict: i --> genotype of mutation at position i (0-based)
    mutation_dict = dict()

    for pos in pos_list:

        orig_allele = ref[pos]  # Original allele
        
        # A list of two alleles: the original and a random other allele
        allele_space = [orig_allele] + random.sample(set(alphabet) - set([orig_allele]), 1)
        
        # Generate mutation genotype
        freq_list = generate_allele_freq(allele_space, 'UNIFORM')
        while True: 
            genotype = [weighted_random(allele_space, freq_list) for i in range(pl)]                        
            # Repeat if every new allele matches the original allele
            if not all([x==orig_allele for x in genotype]): break

        # Mutate position
        for a, s in zip(genotype, seq_list): s[pos] = a

        mutation_dict[pos] = genotype

    seq_list = [''.join(s) for s in seq_list]

    return seq_list, mutation_dict

def generate_normal_genome(ref, mu, pl, distribution='UNIFORM'):
    return mutate_seq(ref, mu, pl, distribution)

def generate_cancer_genome(normal_genome, mu, distribution='UNIFORM'):

    seq_list = []  # List of mutated haploids

    # dict: i --> genotype of mutation at position i (0-based)
    mutation_dict = dict()

    # Independently mutate each haploid of the normal genome
    for h, seq in enumerate(normal_genome):
        s_list, m_dict = mutate_seq(seq, mu, 1, distribution)
        
        seq_list.extend(s_list)

        for pos, genotype in m_dict.items():
            if not mutation_dict.has_key(pos):
                mutation_dict[pos] = [s[pos] for s in normal_genome]
            mutation_dict[pos][h] = genotype[0]

    return seq_list, mutation_dict

def main():
    parser = argparse.ArgumentParser(description='Generate normal or cancer genomes by mutating a reference sequence')
    parser.add_argument('--ref', required=True, help='FASTA file of a contiguous reference sequence')
    parser.add_argument('--k', default=['n1, c1'], nargs='+', required=True, help="List of genome types to generate.  Each genome type is denoted by a two-letter word where the first letter is either 'n' (normal) or 'c' (cancer).  The second letter is a positive integer indicating which individual the genome comes from.  For example, 'n1 c1 n2' means generating three genomes where the first two genomes are a normal and cancer genome from the same individual, and where the third genome is a normal one taken from a different individual.  As another example, 'n1 n2 n3 n4' means generating four normal genomes, each of which comes from a different individual.")
    parser.add_argument('--pl', required=True, type=int, help='Genome ploidy')
    parser.add_argument('--germ-mu', type=float, required=True, help='Germline mutation rate, i.e. the probability that any haploid in a normal genome differs from the reference at a particular site.')
    parser.add_argument('--som-mu', type=float, required=True, help='Somatic x mutation rate, i.e. the probability that any haploid in a cancerous genome differs from the reference at a particular site.')
    parser.add_argument('--output-dir', default=os.getcwd(), help="Directory of output files (Default: current directory)")
    args = parser.parse_args()

    # Check params
    for x in args.k:
        assert x[0] in ['n', 'c'], "Bad format for --k.  Genome type must either be 'n' for normal or 'c' for cancer"        
        try: int(x[1])
        except ValueError: raise Exception("Bad format for --k. Use positive integer to denote individual")
        assert len(x)==2, "Bad format for --k.  Enter two-letter genome type descriptors."
    # if not os.path.isdir(args.output_dir):
    #     query = "%s does not exit.\nPush <Enter> to create it, or Ctrl-C to exit" % (os.path.abspath(args.output_dir))
    #     x = input(query)
        # try:
        #     os.makedirs(args.output_dir)
        # except:
        #     print 'wacka'
        #     raise Exception("Couldn't create %s" % os.path.abspath(args.output_dir))

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # Read reference genomic region
    header = open(args.ref).read().split('\n')[0]
    assert '>' in header, 'First line of FASTA file needs to be the header'
    header = header[1:]  # Remove '>'
    ref_seq = ''.join(open(args.ref).read().split('\n')[1:])  
    individuals = set([x[1] for x in args.k])


    normal_dict = dict()  # dict: i --> normal genome of individual i
    germ_dict   = dict()  # dict: i --> germline mutations in normal genome of individual i
    canc_dict   = dict()  # dict: i --> somatic mutatinos in cancer genome of individual i

    # Sort genome types by indicies of individuals (doesn't affect
    # correctness but allows for easier memory management)
    args.k.sort(key=lambda x: x[1])

    ### Construct genomes ###
    i = None
    for d in args.k:
        
        # Delete previous individual's normal genome which is no longer necessary
        if (i is not None) and d[1]!=i: del normal_dict[i]

        t = d[0]  # Type of genome: normal or cancer
        i = d[1]  # Index of individual        
        
        # Construct normal genome of individual if it's requested or
        # needed to construct cancer genome
        if t=='n' or (t=='c' and not normal_dict.has_key(i)):  
            seq_list, mutation_dict = generate_normal_genome(ref_seq, args.germ_mu, args.pl)
                        
            normal_dict[i] = seq_list
            germ_dict[i]   = mutation_dict
            
            # Write mutations to file
            tab = '\n'.join(['%s\t%s\t%s' % (pos, ref_seq[pos], '\t'.join([x+'\tGERM' for x in genotype])) \
                       for pos, genotype in sorted(mutation_dict.items(), key=lambda x:x[0])])
            f = open(os.path.join(args.output_dir,d+'.var'), 'w')
            f.write('pos(0-based)\tref.allele\t%s\n' % '\t'.join(['hap_%i.allele\thap_%i.type' % (j,j) for j in range(args.pl)]))
            f.write(tab + '\n')
            f.close()

        # Construct cancer genome by mutating normal genome of the same individual        
        elif t=='c':
            seq_list, som_mut_dict = generate_cancer_genome(normal_dict[i], args.som_mu, args.pl)

            ### Write mutations to file ###
            germ_mut_dict = germ_dict[i]            

            t = dict()  # dict: i --> List of the mutation type ('GERM' or 'SOM') of the genotype at position i        
            m = dict()  # dict: i --> Genotype of mutation (with respect to reference) at position i

            for pos in sorted(set(germ_mut_dict.keys()) | set(som_mut_dict.keys())):
                if som_mut_dict.has_key(pos) and germ_mut_dict.has_key(pos):
                    m[pos] = som_mut_dict[pos]
                    t[pos] = ['SOM' if x!=y else ('GERM' if x!=ref_seq[pos] else 'REF') \
                                  for x,y in zip(germ_mut_dict[pos], som_mut_dict[pos])]
                elif som_mut_dict.has_key(pos):
                    m[pos] = som_mut_dict[pos]
                    t[pos] = ['SOM' if y!=ref_seq[pos] else 'REF' for y in som_mut_dict[pos]]
                elif germ_mut_dict.has_key(pos):
                    m[pos] = germ_mut_dict[pos]
                    t[pos] = ['GERM' if x!=ref_seq[pos] else 'REF' for x in germ_mut_dict[pos]]
                else: raise Exception('Yikes!')
            
            tab = '\n'.join(['%s\t%s\t%s' % (pos, ref_seq[pos], '\t'.join(['\t'.join(x) for x in zip(m[pos], t[pos])])) \
                                 for pos in sorted(t.keys())])
            f = open(os.path.join(args.output_dir,d+'.var'), 'w')
            f.write('pos(0-based)\tref.allele\t%s\n' % '\t'.join(['hap_%i.allele\thap_%i.type' % (j,j) for j in range(args.pl)]))
            f.write(tab + '\n')
            f.close()

        # Write genome to FASTA file
        fasta = '\n'.join([format_fasta('%s | haploid %s' % (header,j+1), s) for j, s in enumerate(seq_list)])
        open(os.path.join(args.output_dir, d+'.fa'), 'w').write(fasta)

if __name__=='__main__':
    main()   
