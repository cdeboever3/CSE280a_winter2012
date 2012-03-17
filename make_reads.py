# Chris DeBoever
# cdeboeve@ucsd.edu

# This script generates a pbs or shell script that preprocesses data from Vancouver by trimming adapter sequences, removing low complexity reads, etc.

import sys, argparse, os, pdb
from RNAseq_tools import *

### helper functions ###

def print_script(inN,outD,basename,scriptN,shell,tempD,RNAseq_pipelineD,pbs_threads,jname):
    """Print pbs/shell script to file
    
    Keyword arguments:
    inN -- the input bam file
    outD -- the output directory for the script and script results
    basename -- used for naming directories and files
    shell -- boolean indicating whether we are creaing a shell script
    tempD -- directory where temporary files are stored
    RNAseq_pipelineD -- directory containing scripts that might be used/needed
    pbs_threads -- number of threads to use for PBS script
    jname -- job name used for pbs and for naming files/directories
    """

    ### magic variables ###
    min_len = 20 # minimum length for reads to be kept
    adapters = '-b ' + ' -b '.join(['AAGCAGTGGTATCAACGCAGAGTACGCGGG','AAGCAGTGGTATCAACGCAGAGTAC','AAGCAGTGGTATCAACGCAGAGT']) # adapters to trim with cutadapt
    cutadapt_n = 2 # the maximum number of times to trim adapters for a read
    lc_method = 'dust'
    lc_threshold = '7'
   
    ### file names ###
    # these are files that will be created as this pbs/shell script runs. Some of these are temporary files, some are output files that will be copied to the output directory at the end, and some are ouput files written directly to the output directory (if they are small)
    input_R1                = basename + '_R1.fastq' # input R1 fastq file. The failed filter reads have alredy been removed 
    input_R2                = basename + '_R2.fastq' # input R2 fastq file. The failed filter reads have alredy been removed
    cut_R1_fastq            = basename + '_R1_tag_cleaned.fastq' # R1 fastq file after trimming adapters
    cut_R2_fastq            = basename + '_R2_tag_cleaned.fastq' # R2 fastq file after trimming adapters
    cut_report_R1           = basename + '_R1_cutadapt_report.txt' # R1 cutadapt report
    cut_report_R2           = basename + '_R2_cutadapt_report.txt' # R2 cutadapt report
    prinseq_good_R1_prefix  = basename + '_R1_prinseq_good' # prefix for good R1 reads from prinseq
    prinseq_good_R2_prefix  = basename + '_R2_prinseq_good' # prefix for good R2 reads from prinseq
    prinseq_bad_R1_prefix   = basename + '_R1_prinseq_bad' # prefix for bad R1 reads from prinseq
    prinseq_bad_R2_prefix   = basename + '_R2_prinseq_bad' # prefix for bad R2 reads from prinseq
    prinseq_good_R1         = basename + '_R1_prinseq_good.fastq' # good R1 reads from prinseq
    prinseq_good_R2         = basename + '_R2_prinseq_good.fastq' # good R2 reads from prinseq
    prinseq_bad_R1          = basename + '_R1_prinseq_bad.fastq' # bad R1 reads from prinseq
    prinseq_bad_R2          = basename + '_R2_prinseq_bad.fastq' # bad R2 reads from prinseq
    cleaned_R1              = basename + '_R1_cleaned.fastq' # the final R1 fastq file
    cleaned_R2              = basename + '_R2_cleaned.fastq' # the final R2 fastq file
    orphan_reads            = basename + '_orphans.fastq' # reads whose mate was removed
    cleaned_R1_gz           = basename + '_R1_cleaned.fastq.gz' # the final R1 fastq file gzipped
    cleaned_R2_gz           = basename + '_R2_cleaned.fastq.gz' # the final R2 fastq file gzipped 
    orphan_reads_gz         = basename + '_orphans.fastq.gz' # reads whose mate was removed gzipped
    
    removeL = [input_R1,inputR2,cut_R1_fastq,cut_R2_fastq,prinseq_good_R1,prinseq_good_R2,prinseq_bad_R1,prinseq_bad_R2,cleaned_R1,cleaned_R2,orphan_reads] # these files should be removed at the end of the script. This is important when tempdir==outdir so you don't have a bunch of extra files lying around

    ### some other useful stuff ###
    current_time = get_timestamp() # timestamp for script creation
    git_info = get_git_info('{0}/.git'.format(RNAseq_pipelineD)) # RNAseq_pipelineD is the directory with the RNAseq pipeline scripts
   
    ### let's write the script ###
    scriptF = open(scriptN,'w')
    
    print >>scriptF, '#!/bin/bash'
    print >>scriptF, ''
    print >>scriptF, '# This script preprocesses Vancouver RNA-seq data'
    print >>scriptF, ''
    print >>scriptF, tracking_info(current_time,' '.join(sys.argv[0:]),RNAseq_pipelineD,git_info)
    print >>scriptF, ''
    if pbs:
        print >>scriptF, '#PBS -N ' + jname
        print >>scriptF, '#PBS -l nodes=1:ppn={0}'.format(pbs_threads)
        print >>scriptF, '#PBS -o {0}/{1}.out'.format(outD,jname)
        print >>scriptF, '#PBS -e {0}/{1}.err'.format(outD,jname)
        print >>scriptF, ''
        print >>scriptF, 'rm -r ' + tempD # remove the tempD if it already exists to avoid problems
        print >>scriptF, ''
    print >>scriptF, 'mkdir -p ' + tempD
    print >>scriptF, 'mkdir -p ' + outD
    print >>scriptF, ''
    print >>scriptF, 'cd ' + tempD
    print >>scriptF, ''
    # copy needed scripts to this directory to avoid version discrepancies if the scripts are updated
    print >>scriptF, 'rsync -avz {0}/remove_orphan_reads.py .'.format(RNAseq_pipelineD)
    print >>scriptF, ''
    # copy bam file to this directory
    print >>scriptF, 'rsync -avz {0} .'.format(inN)
    print >>scriptF, ''
    # convert bam file to fastq files
    print >>scriptF, 'java -Xmx2g -Djava.io.tmpdir={0} -jar /raid/software/src/picard-tools-{1}/SamToFastq.jar VALIDATION_STRINGENCY=SILENT I={2} FASTQ={3} SECOND_END_FASTQ={4} '.format(tempD,picard_version,input_R1,input_R2)
    print >>scriptF, ''
    # collect fastqc statistics before cleaning
    print >>scriptF, '/raid/software/src/fastqc_v{2}/fastqc -o {0} --extract {1} &'.format(outD,input_R1,fastqc_version)
    print >>scriptF, '/raid/software/src/fastqc_v{2}/fastqc -o {0} --extract {1} &'.format(outD,input_R2,fastqc_version)
    print >>scriptF, ''
    # remove adapters with cutadapt
    print >>scriptF, '/raid/software/src/cutadapt-{2}/cutadapt -m {4} -n {6} {5} {0} -o {1} > {4}/{3} &'.format(input_R1,cut_R1_fastq,cut_report_R1,outD,min_len,adapters,cutadapt_n)
    print >>scriptF, '/raid/software/src/cutadapt-{2}/cutadapt -m {4} -n {6} {5} {0} -o {1} > {4}/{3}'.format(input_R2,cut_R2_fastq,cut_report_R2,outD,min_len,adapters,cutadapt_n)
    print >>scriptF, ''
    print >>scriptF, 'wait'
    print >>scriptF, ''
    # remove low complexity and short reads with prinseq
    print >>scriptF, '/raid/software/src/prinseq-lite-{0}/prinseq-lite -verbose -lc_method {1} -lc_threshold {2} -min_len {3} -fastq {4} -out_good {5} -out_bad {6} &'.format(prinseq_lite_version,lc_method,lc_threshold,min_len,cut_report_R1,prinseq_good_R1_prefix,prinseq_bad_R1_prefix)
    print >>scriptF, '/raid/software/src/prinseq-lite-{0}/prinseq-lite -verbose -lc_method {1} -lc_threshold {2} -min_len {3} -fastq {4} -out_good {5} -out_bad {6} &'.format(prinseq_lite_version,lc_method,lc_threshold,min_len,cut_report_R2,prinseq_good_R2_prefix,prinseq_bad_R2_prefix)
    print >>scriptF, ''
    print >>scriptF, 'wait'
    print >>scriptF, ''
    print >>scriptF, 'python remove_orphan_reads.py {} {} {} {} {}'.format(prinseq_good_R1,prinseq_good_R2,cleaned_R1,cleaned_R1,orphan_reads)
    print >>scriptF, ''
    print >>scriptF, '/raid/software/src/fastqc_v{2}/fastqc -o {0} --extract {1} &'.format(outD,cleaned_R1,fastqc_version)
    print >>scriptF, '/raid/software/src/fastqc_v{2}/fastqc -o {0} --extract {1} &'.format(outD,cleaned_R2,fastqc_version)
    print >>scriptF, ''
    print >>scriptF, 'wait'
    print >>scriptF, ''
    print >>scriptF, 'gzip {0}_R1_cleaned.fastq &'.format(basename)
    print >>scriptF, 'gzip {0}_R2_cleaned.fastq &'.format(basename)
    print >>scriptF, 'gzip {0}_orphans.fastq'.format(basename)
    print >>scriptF, 'gzip {0} &'.format(cleaned_R1)
    print >>scriptF, 'gzip {0} &'.format(cleaned_R2)
    print >>scriptF, 'gzip {0}'.format(orphan_reads)
    print >>scriptF, ''
    print >>scriptF, 'rm {0}'.format(inN.split('/')[-1])
    if tempD != outD: # if the tempD isn't the outD, we need to copy files to the outD
        print >>scriptF, ''
        print >>scriptF, 'wait'
        print >>scriptF, ''
        print >>scriptF, 'rsync -av *.py {0}'.format(outD) # copy over scripts for record keeping purposes
        # copy final fastq files to output directory
        print >>scriptF, 'rsync -avz {1}_R1_cleaned.fastq.gz {0} &'.format(outD,basename) 
        print >>scriptF, 'rsync -avz {1}_R2_cleaned.fastq.gz {0} &'.format(outD,basename)
        print >>scriptF, 'rsync -avz {1}_orphans.fastq.gz {0}'.format(outD,basename)
        print >>scriptF, ''
        print >>scriptF, 'wait'
        print >>scriptF, ''
        print >>scriptF, 'rm -r {0}'.format(tempD)
        print >>scriptF, ''
    if tempD == outD: # if the tempD is the outD, we want to remove some temporary files
        print >>scriptF, 'rm {}'.format(' '.join(removeL))
        print >>scriptF, ''
    print >>scriptF, ''
    scriptF.close()

def main():
    ### magic variables ###
    jname_suffix = '_reads' # used for naming jobs and files
    pbs_threads = 8 # 6 threads here is more than enough because nothing is multi-threaded. However, we want to take up some space on the node just so we don't have too much IO on the hdd that slows everything down
    # the following variables can be set at the command line
    # if a variable is set to 'default', then it relies on another variable that is set required at the command line
    fastaL      = ['/raid/development/cdeboever/CSE280a/n1.fa', '/raid/development/cdeboever/CSE280a/c1.fa']
    covL        = [30]
    alphaL      = [0.05]
    outD        = 'default'
    basename    = 'default'
    scriptN     = 'default'
    tempD       = 'default'

    ### gather arguments from command line ###
    parser = argparse.ArgumentParser(description='This script makes a PBS or shell for making reads. The PBS script is saved in outdir.')
    parser.add_argument('-f', metavar='fasta_list', default=fastaL, nargs='+', help='Reference fasta files. Default: {0}'.format(' '.join(fastaL)))
    parser.add_argument('-c', metavar='coverage_list', default=covL, nargs='*', help='Desired coverage levels. Default: {0}'.format(' '.join(covL)))
    parser.add_argument('-a', metavar='alpha_list', default=alphaL, nargs='*', help='Desired alpha values. Default: {0}'.format(' '.join(alphaL)))
    parser.add_argument('-o', metavar='outdir', default=outD, help='The directory where the output directory should be stored. If not specified, the output folder will be made in the current directory. The directory will be made if it does not exist. Default: [basename]_{0}'.format(jname_suffix))
    parser.add_argument('-b', metavar='basename', default=basename, help='Basename for naming output files. The basename will be read from the path of the input bam file if not specified.')
    parser.add_argument('-s', metavar='script_name', default=scriptN, help='Output pbs script (default: [basename]_{0}.[pbs][sh])'.format(jname_suffix))
    parser.add_argument('-t', metavar='tempdir', default=tempD, help='Directory for storing temporary files. This is useful when making a shell script but likely should not be specified when making PBS scripts. WARNING: when this option is used, this tempdir is deleted at the end of the shell script.')
    parser.add_argument('--shell', action='store_true', help='Output a shell script. WARNING: temporary files will be stored in the output directory. Choose the output directory wisely or specify the temp directory.')
    parser.add_argument('--debug', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
   
    fastaL      = args.f        # list of fasta files to be used
    covL        = args.c        # list of coverages to be used
    alphaL      = args.a        # list of alphas to be used
    outD        = args.o        # output directory
    basename    = args.b        # basename used for naming output directory and files
    scriptN     = args.s        # name of output PBS/shell script
    shell       = args.shell    # boolean indicating whether to make shell script as opposed to PBS script
    tempD       = args.t        # directory where temporary files should be stored
    debug       = args.debug
    
    ### let's set some variables that we'll use throughout ###
    
    # get the git directory with the code
    RNAseq_pipelineD = os.path.dirname(os.path.realpath(__file__))
    
    # get the full path of the read directory
    inN = os.path.realpath(inN)
    
    # a couple booleans for convenience
    if shell: # true if we are making a shell script
        pbs = False
    else:
        pbs = True # true if we are making a pbs script
        shell = False
    
    # set the basename
    if basename == 'default': # used for naming output files, folders, job names, etc.
        basename = os.path.dirname(inN).split('/')[-1]
    
    # set the job name
    jname = basename + jname_suffix # job name for PBS job. We'll use this for naming file and folders as well
    
    # set the output directory
    if outD == 'default':
        outD = '{0}/{1}'.format(os.getcwd(),jname)
    else:
        outD = '{0}/{1}'.format(os.path.realpath(outD),jname)
    
    # set the temp directory
    if tempD == 'default': # true if we need to set the temp directory
        if shell:
            tempD = outD # if you are writing a shell script, we'll assume the output directory is where you want temporary files to be stored
        else:
            tempD = '/scratch/' + jname # for pbs scripts, temp files are stored on the /scratch disk
    
    # set the pbs/shell script name
    if scriptN == 'default':
        if pbs:
            scriptN = '{0}/{1}.pbs'.format(outD,jname)
        else:
            scriptN = '{0}/{1}.sh'.format(outD,jname)
    
    # make the output directory
    try:
        os.makedirs(outD)
    except OSError:
        pass
    
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
   
    print_script(inN,outD,basename,scriptN,shell,tempD,RNAseq_pipelineD,pbs_threads,jname)

    sys.exit(0)

if __name__ == '__main__':
    main()
