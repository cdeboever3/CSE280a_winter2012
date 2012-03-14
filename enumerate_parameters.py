import sys, math, os
from partition_list import partition_list
import subprocess, multiprocessing

def run(cmd):
    '''Spawn a shell command and wait for it to finish.  Raise an
    exception if the command was not successful, as indicated by a
    non-zero return code.'''
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    try:
        assert p.returncode==0
    except:
        stdout, stderr = p.communicate()
        print 'Command failure:', cmd
        print 'Command stdout:', stdout
        print 'Command stderr:', stderr
    
def make_script(param1_list, param2_list, param1_type, param2_type, folder, cores):
    '''
    Run generate_sample_tsv_complex.tsv and
    mixed_variant_calling.py for all pairs (param1, param2) from
    <param1_list> and <param2_list>.  Output files to
    <folder>. Distribute processing across <cores> processes.
    
    @param1_type : Name of parameter listed in param1_list
    @param2_type : Name of parameter listed in param2_list
    
    Currently only supports param1
    '''

    cmd_list = []  # A list of shell commands to run

    for param1 in param1_list:
        for param2 in param2_list:

            # Create command
            if param1_type=='coverage' and param2_type=='alpha':
                c, a = param1, param2
                tN = 'c_{0}_a_{1}.tsv'.format(c,a)
                rN = 'c_{0}_a_{1}_truth.txt'.format(c,a)
                gN = 'c_{0}_a_{1}_genotypes.txt'.format(c,a)
                tN, rN, gN = [os.path.join(folder,x) for x in [tN, rN, gN]]
                cmd1 = 'python generate_sample_tsv_complex.tsv {0} -o {1} -ac {2} -sl 1000000 -sL {3} {4}'.format(rN,tN,c,a,1-a)
                cmd2 = 'python mixed_variant_calling.py {2} -er {0} > {1}'.format(e,gN,tN)            
            elif param1_type=='alpha' and param2_type=='error':
                a, e = param1, param2
                tN = 'a_{0}_e_{1}.tsv'.format(a,e)
                rN = 'a_{0}_e_{1}_truth.txt'.format(a,e)
                gN = 'a_{0}_e_{1}_genotypes.txt'.format(a,e)
                tN, rN, gN = [os.path.join(folder,x) for x in [tN, rN, gN]]
                cmd1 = 'python generate_sample_tsv_complex.tsv {0} -o {1} -sl 1000000 -er {2} -sL {3} {4}'.format(rN,tN,e,a,1-a)
                cmd2 = 'python mixed_variant_calling.py {1} -er 0.001 > {0}'.format(gN,tN)
            else:
                raise Exception('Param types %s and %s not yet supported' % (param1_type, param2_type))

            cmd = cmd1 + ' ; ' + cmd2  # Concatenate commands

            cmd_list.append(cmd)

    # Distribute shell commands across multiple processors
    pool = multiprocessing.Pool(processes=cores)
    pool.map(run, cmd_list)

if __name__=='__main__':

    ### Set one of these to be True and the other False
    alpha_vs_error = True
    coverage_vs_alpha = False

    assert alpha_vs_error + coverage_vs_alpha == 1

    cores = 24

    if alpha_vs_error:

        errorL = [0.001, 0.01, 0.05]
        alpha_step = 0.1
        alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]

        folder = 'alpha_vs_error'
        if not os.path.isdir(folder): os.makedirs(folder)

        make_script(alphaL, errorL, 'alpha', 'error', folder, cores)

    elif coverage_vs_alpha:
        covL = [2**i for i in range(2, 8)]
        alpha_step = 0.05
        alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]

        folder = 'coverage_vs_alpha'
        if not os.path.isdir(folder): os.makedirs(folder)

        make_script(covL, alphaL, 'coverage', 'alpha', folder, cores)
