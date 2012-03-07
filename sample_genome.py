#!/usr/bin/python

import os, argparse, subprocess, glob, random, itertools, bisect

def get_fasta_lengths(genome_dir, precomputed_file=None):
    """
    Returns a dictionary of d[chromosome] := chromosome length.  FASTA
    files of chromosomes are read from <genome_dir>.  If
    <precomputed_file> is specified, then (chromosome, chromosome
    length) pairs are read from the file if it exists, otherwise the
    pairs are written to the file in a 2-column tab-delimited format.
    """
    if os.path.isfile(precomputed_file):
        # Read lengths from precomputed file
        d = dict([(x.split()[0], int(x.split()[1])) \
                      for x in open(precomputed_file).read().split('\n') if x!=''])        
    else:
        # Calculate lengths
        fasta_files = glob.glob(os.path.join(genome_dir, '*.fa'))
        d = dict([(os.path.basename(f), len(''.join(open(f).read().split('\n')[1:]))) \
                      for f in fasta_files])

        if precomputed_file is not None:
            # Write lengths to file
            open(precomputed_file, 'w').write( \
                '\n'.join(['\t'.join(map(str, x)) for x in d.items()]) + '\n')
    return d

def weighted_random(a, b):
    """Randomly selects an element a[i] with probability proportional to b[i].
    See  http://docs.python.org/dev/py3k/library/random.html#examples-and-recipes."""
    b = b[:]  # Copy list
    for i in range(len(b))[1:]: b[i] += b[i-1]
    return a[bisect.bisect_right(b, random.random() * b[-1])]

def format_fasta(header, seq, width=50):
    """Prepares a header and a sequence String into FASTA format.  For
    prettier output, a newline is inserted after every <width>
    characters in <seq>."""
    
    # Insert newlines
    seq_lines, i = [], 0
    while (width*i < len(seq)):
        seq_lines.append(seq[width*i: width*(i+1)])
        i += 1
    seq = '\n'.join(seq_lines)
    return '>%s\n%s\n' % (header,seq)

def main():
    parser = argparse.ArgumentParser(description='Randomly selects a contiguous region from a genome.  Removes soft masking by making all characters uppercase.')
    parser.add_argument('--genome-dir', required=True, help='Directory of *.fa files of a genome.')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('--sample-length', default=int(1e8), type=float, help='Length of region to sample (Default: 100 Mb).')
    args = parser.parse_args()

    # Check params
    assert os.path.isdir(args.genome_dir), 'Error: %s is not a valid directory' % args.genome_dir
    try:  args.sample_length = int(args.sample_length)
    except ValueError: raise Exception('Input an integer for --sample-length')

    # Dict: chromosome --> chromosome length
#    length_dict = get_fasta_lengths(args.genome_dir)
    length_dict = get_fasta_lengths(args.genome_dir, precomputed_file='chrom_lengths.txt')

    # Filter for chromosomes with > 100Mb
    length_dict = dict([(k,v) for k,v in length_dict.items() if v > args.sample_length])
    
    # Randomly pick a chromosome with prob proportional to its length
    chrom_list, length_list = map(list, zip(*(length_dict.items())))
    chrom = weighted_random(chrom_list, length_list)

    # Randomly select 100 Mb sequence from chromosomee
    start = random.randint(1, length_dict[chrom] - args.sample_length)
    chrom_seq = ''.join(open(os.path.join(args.genome_dir, chrom)).read().split('\n')[1:])
    sample_seq = chrom_seq[start - 1: start + args.sample_length - 1]
    
    sample_seq = sample_seq.upper()  # Upper case

    # Write sequence to FASTA file
    fasta = format_fasta('%s:%s-%s' % (chrom, start, start+args.sample_length), sample_seq)
    open(args.output, 'w').write(fasta)
        
if __name__=='__main__':
    main()
