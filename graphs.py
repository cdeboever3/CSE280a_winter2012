"""
Tabulates the variant calling results produced by
enumerate_parameters.py.  The tables can then be graphed in R with
graphs.R
"""

import os, glob, argparse

def write_table_header(fh, param_list):
    "Write table header to open file <fh>"
    fh.write('\t' + '\t'.join(map(str,param_list)) + '\n')

def read_genotypes(genotypes_file):
    """
    Input: '_genotypes.txt' or '_truth.txt' file.
    Output: alpha and dictionaries mapping genomic position to cancer and normal genotype
    """

    lines = [x for x in open(genotypes_file).read().split('\n') if x!='']
    alpha = [x.split('\t')[1:3] for x in lines if x[0]=='#' and 'population frequencies' in x]
    assert len(alpha)==1
    alpha = alpha[0]
    lines = [x.split('\t') for x in lines if x[0]!='#']

    # dict: genomic position --> genotype of cancer
    cancer_dict = dict([(x[0]+'.'+x[1], set(x[3].upper())) for x in lines])

    # dict: genomic position --> genotype of normal
    normal_dict = dict([(x[0]+'.'+x[1], set(x[2].upper())) for x in lines])

    return alpha, cancer_dict, normal_dict

def main():

    ### Setup arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input-folder', required=True)
    parser.add_argument('--output-folder', required=True)
    args = parser.parse_args()

    input_folder = args.input_folder
    output_folder = args.output_folder

    # Infer list of parameters
    files = [x for x in glob.glob(os.path.join(input_folder, 'c_*_a_*_e_*_genotypes.txt')) if os.path.getsize(x)>0]
    covL = sorted(list(set([int(x.split('c_')[1].split('_')[0]) for x in files])))
    alphaL = sorted(list(set([float(x.split('c_')[1].split('_')[2]) for x in files])))
    errorL = sorted(list(set([float(x.split('c_')[1].split('_')[4]) for x in files])))

    # Hack: remove 0.05 from errorL
    if 0.05 in errorL: errorL.remove(0.05)

    print covL
    print alphaL
    print errorL

    c_vs_a = False
    
    if c_vs_a:
        param1_list = covL
        param2_list = alphaL
        e = errorL[0]
    else:
        param1_list = alphaL
        param2_list = errorL
        c = covL[0]

    ## Output filenames
    base = os.path.join(output_folder, 'alpha_vs_error')
    tp_file = open(base+'.tp.txt', 'w')        # True positives
    fp_file = open(base+'.fp.txt', 'w')        # False positives
    fn_file = open(base+'.fn.txt', 'w')        # False negatives
    tn_file = open(base+'.tn.txt', 'w')        # True negatives
    alpha_file = open(base+'.alpha.txt', 'w')  # Alpha predictions
    tv_file = open(base+'.tv.txt', 'w')        
    tnv_file = open(base+'.tnv.txt', 'w')
    pv_file = open(base+'.pv.txt', 'w')
    pnv_file = open(base+'.pnv.txt', 'w')

    ## Setup output files
    for x in [tp_file, fp_file, fn_file, tn_file, alpha_file, tv_file, tnv_file, pv_file, pnv_file]:
        write_table_header(x, param2_list)

    ### Enumerate over parameters
    for param1 in param1_list:

        tp_list, fp_list, fn_list, tn_list, alpha_list = [], [], [], [], []

        tv_list, tnv_list, pv_list, pnv_list = [], [], [], []

        for param2 in param2_list:

            if c_vs_a:
                c, a = param1, param2
            else:
                a, e = param1, param2

            # Input filenames
            rN = 'c_{0}_a_{1}_e_{2}_truth.txt'.format(c,a,e)
            gN = 'c_{0}_a_{1}_e_{2}_genotypes.txt'.format(c,a,e)
            rN, gN = [os.path.join(input_folder,x) for x in [rN, gN]]

            # Read files
            true_alpha, true_cancer_dict, true_normal_dict = read_genotypes(rN)
            predict_alpha, predict_cancer_dict, predict_normal_dict = read_genotypes(gN)

            # False and True positives
            p,q,r,s = true_cancer_dict, true_normal_dict, predict_cancer_dict, predict_normal_dict
            fp = len([k for k in r.keys() if r[k]!=s[k] and p[k]==q[k]])
            tp = len([k for k in r.keys() if r[k]!=s[k] and p[k]!=q[k]])
            fn = len([k for k in r.keys() if r[k]==s[k] and p[k]!=q[k]])
            tn = len([k for k in r.keys() if r[k]==s[k] and p[k]==q[k]])

            # # False and True negatives
            # p,q,r,s = predict_cancer_dict, predict_normal_dict, true_cancer_dict, true_normal_dict

            # actual alpha - predicted alpha
            alpha_error = float(true_alpha[0]) - float(predict_alpha[0])

            true_var = [k for k in true_cancer_dict.keys() if true_cancer_dict[k]!=true_normal_dict[k]]
            true_non_var = [k for k in true_cancer_dict.keys() if true_cancer_dict[k]==true_normal_dict[k]]
            assert len(true_var) + len(true_non_var) == len(true_cancer_dict.keys())

            predict_var = [k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]!=predict_normal_dict[k]]
            predict_non_var = [k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]==predict_normal_dict[k]]
            assert len(predict_var) + len(predict_non_var) == len(predict_cancer_dict.keys())

            assert tp + fp == len(predict_var)
            assert tn + fn == len(predict_non_var)

            alpha_list.append(alpha_error)
            tp_list.append(tp)
            fp_list.append(fp)
            fn_list.append(fn)
            tn_list.append(tn)
            tv_list.append(len(true_var))
            tnv_list.append(len(true_non_var))
            pv_list.append(len(predict_var))
            pnv_list.append(len(predict_non_var))

        print >>tp_file, str(param1) + '\t' + '\t'.join(map(str,tp_list))
        print >>fp_file, str(param1) + '\t' + '\t'.join(map(str,fp_list))
        print >>fn_file, str(param1) + '\t' + '\t'.join(map(str,fn_list))
        print >>tn_file, str(param1) + '\t' + '\t'.join(map(str,tn_list))
        print >>alpha_file, str(param1) + '\t' + '\t'.join(map(str,alpha_list))
        print >>tv_file, str(param1) + '\t' + '\t'.join(map(str,tv_list))
        print >>tnv_file, str(param1) + '\t' + '\t'.join(map(str,tnv_list))
        print >>pv_file, str(param1) + '\t' + '\t'.join(map(str,pv_list))
        print >>pnv_file, str(param1) + '\t' + '\t'.join(map(str,pnv_list))

    tp_file.close()
    fp_file.close()
    fn_file.close()
    tn_file.close()
    alpha_file.close()
    tv_file.close()
    tnv_file.close()
    pv_file.close()
    pnv_file.close()

if __name__=='__main__':
    main()
