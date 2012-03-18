import sys, math, os
import subprocess, multiprocessing, argparse

def run_cmd(cmd, debug):
    '''Spawn a shell script <cmd> and wait for it to finish.  Raise an
    exception if the command was not successful, as indicated by a
    non-zero return code.'''

    if debug: print >>sys.stderr, 'Starting command:', cmd
    
    # Run and wait for shell script
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    try:        
        assert p.returncode==0
    except:
        # Script was not successful
        if debug:
            stdout, stderr = p.communicate()
            print >>sys.stderr, 'Command failure:', cmd
            print >>sys.stderr, 'Command stdout:', stdout
            print >>sys.stderr, 'Command stderr:', stderr
    if debug: print 'Finished command:', cmd

def run_cmd_star(x):
    return run_cmd(*x)

def make_script(param_list, folder, cores, debug=False):
    """
    Run generate_sample_tsv_complex.tsv and mixed_variant_calling.py
    for all parameters in <param_list>.  Output files to
    <folder>. Distribute processing across <cores> processes.
    
    @param_list : A list of 3-tuples (alpha, error, coverage).
    """

    cmd_list = []  # A list of shell commands to run

    for a, e, c in param_list:

        # Output filenames
        tN = 'c_{0}_a_{1}_e_{2}.tsv'.format(c,a,e)
        rN = 'c_{0}_a_{1}_e_{2}_truth.txt'.format(c,a,e)
        gN = 'c_{0}_a_{1}_e_{2}_genotypes.txt'.format(c,a,e)
        tN, rN, gN = [os.path.join(folder,x) for x in [tN, rN, gN]]
            
        # Create command
        cmd1 = 'python generate_sample_tsv_complex.tsv {0} -o {1} -ac {2} -sl 1000000 -sL {3} {4} -er {5}'.format(rN,tN,c,a,1-a,e)
        cmd2 = 'python mixed_variant_calling.py {2} -er {0} > {1}'.format(e,gN,tN)                        

        cmd = cmd1 + ' ; ' + cmd2  # Concatenate commands

        cmd_list.append(cmd)

    # Distribute commands across multiple processors
    pool = multiprocessing.Pool(processes=cores)
    #for x in [(cmd, debug) for cmd in cmd_list][:10]: print x
    pool.map(run_cmd_star, [(cmd, debug) for cmd in cmd_list])

if __name__=='__main__':

    # Determines which parameters to vary
    switch = False

    if switch:
        ### Vary alpha and error rate
        errorL = [0.001, 0.01, 0.02, 0.05]
        alpha_step = 0.1
        alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]
        covL = [50]

        folder = 'alpha_vs_error'
    else:
        ### Vary coverage and alpha
        errorL = [.001]
        covL = [2**i for i in range(2, 8)]
        alpha_step = 0.05
        alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]

        folder = 'coverage_vs_alpha'

    cores = 24

    if not os.path.isdir(folder): os.makedirs(folder)

    param_list = [(a,e,c) for a in alphaL for e in errorL for c in covL]

    make_script(param_list, folder, cores)
