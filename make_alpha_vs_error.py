# Chris DeBoever
# cdeboeve@ucsd.edu

import sys, argparse, os, pdb
from RNAseq_tools import *

parser = argparse.ArgumentParser(description='This script makes a pbs or shell script for simulating reads and running our variant calling pipeline.')
parser.add_argument('genome_fastas', nargs='+', help='A list of genome fasta files.')
parser.add_argument('-pL', default=ploidyL, type=int, nargs='+', help='A list of ploidies. Each entry in the list represents the anticipated ploidy of a subpopulation. For instance, if you expect two diploid subpopulations and one triploid subpopulation, enter 2 2 3. Default: {0}'.format(' '.join([str(x) for x in ploidyL])))
parser.add_argument('-outdir', default='default', help='The directory where the output directory should be stored. The directory will be made if it does not exist.')
parser.add_argument('-basename', default = 'default', help='Basename for naming output files.')
parser.add_argument('-o', metavar='script_name', default='default', help='Output script.')
parser.add_argument('-shell', action='store_true', help='Output a shell script. WARNING: temporary files will be stored in the output directory. Choose the output directory wisely or specify the temp directory.')
parser.add_argument('-tempdir', default='default', help='Directory for storing temporary files. This is useful when making a shell script but likely should not be used when making PBS scripts. WARNING: this tempdir is deleted at the end of the shell script.')
parser.add_argument('-d', action='store_true', help='Enable python debugger.')

args = parser.parse_args()

alignD      = args.aligndir
outD        = args.outdir
basename    = args.basename
scriptN     = args.o
shell       = args.shell
tempD       = args.tempdir
debug       = args.d

### let's set some variables that we'll use throughout ###

# get the git directory with the code
RNAseq_pipelineD = os.path.dirname(os.path.realpath(__file__))

# get the full path of the read directory
alignD = os.path.realpath(alignD)

# a couple booleans for convenience
if shell: # true if we are making a shell script
    pbs = False
else:
    pbs = True # true if we are making a pbs script
    
# remove trailing / from alignD if present
if alignD[-1] == '/':
    alignD = alignD[:-1]

# set the basename
if basename == 'default': # used for naming output files, folders, job names, etc.
    basename = alignD.split('/')[-2]

# set the job name
jname = basename + '_gsnap_process' # job name for PBS job. We'll use this for naming file and folders as well

# set the output directory
if outD == 'default':
    outD = '{0}/{1}'.format(os.path.dirname(alignD),jname)
else:
    outD = '{0}/{1}'.format(os.path.realpath(outD),jname)

# set the temp directory
if tempD == 'default': # true if we need to set the temp directory
    if shell:
        tempD = outD # if you are writing a shell script, we'll assume the output directory is where you want temporary files to be stored
    else:
        tempD = '/scratch/' + jname # for pbs scripts, temp files are stored on the /scratch disk

# set the pbs/shell script name
if pbs:
    scriptN = '{0}/{1}.pbs'.format(outD,jname)
else:
    scriptN = '{0}/{1}.sh'.format(outD,jname)

# make the output directory
try:
    os.makedirs(outD)
except OSError:
    pass

metricD = 'metric_files' # name of folder where metric files will be stored in output directory
codeD = 'code' # name of folder where code will be stored in output directory

### time to write the pbs/shell script ###

# check to see if this pbs/shell script already exists so we don't accidentally overwrite
status = ''
if os.path.isfile(scriptN):
    while status != 'c' or 'q':
        status = raw_input('PBS/shell script already exists, enter c to continue or q to quit.\n')
        if status == 'c':
            break
        elif status == 'q':
            sys.exit(0)

current_time = get_timestamp()
git_info = get_git_info('{0}/.git'.format(RNAseq_pipelineD)) # RNAseq_pipelineD is the directory with the RNAseq pipeline scripts

scriptF = open(scriptN,'w')

print >>scriptF, '#!/bin/bash'
print >>scriptF, ''
print >>scriptF, '# This script processes a GSNAP alignment of RNA-seq data'
print >>scriptF, ''
print >>scriptF, tracking_info(current_time,' '.join(sys.argv[0:]), RNAseq_pipelineD, git_info)
print >>scriptF, ''
if pbs:
    print >>scriptF, '#PBS -N ' + jname
    print >>scriptF, '#PBS -l nodes=1:ppn=8'
    print >>scriptF, '#PBS -o ' + outD + '/' + jname + '.out'
    print >>scriptF, '#PBS -e ' + outD + '/' + jname + '.err'
    print >>scriptF, ''
    print >>scriptF, 'rm -r ' + tempD # remove the tempD if it already exists to avoid problems
print >>scriptF, 'mkdir -p ' + tempD
print >>scriptF, 'mkdir -p ' + outD
print >>scriptF, ''
print >>scriptF, 'cd ' + tempD
print >>scriptF, ''
# copy needed scripts to this directory to avoid version discrepancies if the scripts are updated
print >>scriptF, 'rsync -avz {0}/RNAseq_tools.py .'.format(RNAseq_pipelineD)
print >>scriptF, 'rsync -avz {0}/remove_chimeric_alignments.py .'.format(RNAseq_pipelineD)
print >>scriptF, 'rsync -avz {0}/mark_mult_duplicates.py .'.format(RNAseq_pipelineD)
print >>scriptF, 'rsync -avz {0}/gsnap_flagstat.py .'.format(RNAseq_pipelineD)
print >>scriptF, 'rsync -avz {0}/count_reads.py .'.format(RNAseq_pipelineD)
print >>scriptF, ''
# unique alignments
uniqL = ['.concordant_uniq', \
'.paired_uniq_inv', \
'.paired_uniq_long', \
'.paired_uniq_scr', \
'.unpaired_uniq']
for fN in halfL:
    print >>scriptF, 'rsync -avz {0}/{1}{2}.bam . &'.format(alignD,basename,fN)
print >>scriptF, ''
print >>scriptF, 'wait'
print >>scriptF, ''
for fN in halfL:
    print >>scriptF, 'python remove_chimeric_alignments.py {0}{1}.bam -o {0}{1}_no_chim.bam -RF {0}{1}_chim.bam &'.format(basename,fN,RNAseq_pipelineD)
print >>scriptF, ''
print >>scriptF, 'wait'
print >>scriptF, ''
# for fN in halfL:
#     print >>scriptF, 'rsync -avz {0}{1}_no_chim.bam {0}{1}_chim.bam {2} &'.format(basename,fN,outD)

# translocation alignments
translocL = ['.concordant_transloc', \
'.halfmapping_transloc', \
'.unpaired_transloc']
print >>scriptF, ''
for fN in translocL:
    print >>scriptF, 'rsync -avz {0}/{1}{2}.bam . &'.format(alignD,basename,fN)
print >>scriptF, 'rsync -avz {0}/{1}{2}.bam . &'.format(alignD,basename,'.nomapping')
print >>scriptF, ''
print >>scriptF, 'wait'
print >>scriptF, ''
print >>scriptF, 'python gsnap_flagstat.py {0} -basename {2} > {1}/{4}/{2}_gsnap_flagstat.dat'.format(tempD,outD,basename,RNAseq_pipelineD,metricD) 
print >>scriptF, ''
if pbs or tempD != outD:
    print >>scriptF, ''
    print >>scriptF, 'rsync -avz *py {0}/{1}'.format(outD,codeD)
    print >>scriptF, ''
    print >>scriptF, 'wait'
    print >>scriptF, ''
    print >>scriptF, 'rm -r ' + tempD
print >>scriptF, ''
scriptF.close()

if __name__ == '__main__':
    main()
