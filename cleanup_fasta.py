# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, pdb, random, re
import numpy as np

### helper functions ###

rn = re.compile('N')

def main():

    parser = argparse.ArgumentParser(description='This script takes a fasta file and makes all bases upper cases and replaces N\'s with random nucleotides.')
    parser.add_argument('infile', help='Input fasta file.')
    parser.add_argument('-o', nargs='?', type=argparse.FileType('w'),default=sys.stdout, help='Output file. Default: standard out')
    parser.add_argument('-d', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
  
    inN     = args.infile
    outF    = args.o
    debug   = args.d

    inF = open(inN,'r')
    line = inF.readline().strip()
    while line != '':
        if line[0] != '>':
            line = line.upper()
            while 'N' in line:
                line = rn.sub(random.choice(['A','C','G','T']),line,1)

        print >>outF, line
        line = inF.readline().strip()

if __name__=='__main__':
    main()   
