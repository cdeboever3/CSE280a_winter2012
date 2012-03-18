"""Tabulates the variant calling results produced by
enumerate_parameters.py.  The tables can then be graphed in R with
graphs.R"""

import os, glob
    
# List of parameters
covL = [2**i for i in range(2, 8)]
alpha_step = 0.05
alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]

param1_list = covL
param2_list = alphaL

# Output filenames
folder = 'coverage_vs_alpha'
base = os.path.join('Figures', folder)
tp_file = open(base+'.tp.txt', 'w')        # True positives
fp_file = open(base+'.fp.txt', 'w')        # False positives
fn_file = open(base+'.fn.txt', 'w')        # False negatives
tn_file = open(base+'.tn.txt', 'w')        # True negatives
alpha_file = open(base+'.alpha.txt', 'w')  # Alpha predictions

tv_file = open(base+'.tv.txt', 'w')        
tnv_file = open(base+'.tnv.txt', 'w')
pv_file = open(base+'.pv.txt', 'w')
pnv_file = open(base+'.pnv.txt', 'w')

def write_table_header(fh):
    "Write table header"
    fh.write('\t' + '\t'.join(map(str,param2_list)) + '\n')

for x in [tp_file, fp_file, fn_file, tn_file, alpha_file, tv_file, tnv_file, pv_file, pnv_file]:
    write_table_header(x)

### Enumerate over parameters
fp_list, fn_list, alpha_list = [], [], []

#for a, c, e in param_list:

e = .001

# param1_list = [16]
# param2_list = [0.7]

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

for param1 in param1_list:

    tp_list, fp_list, fn_list, tn_list, alpha_list = [], [], [], [], []

    tv_list, tnv_list, pv_list, pnv_list = [], [], [], []

    for param2 in param2_list:

        c, a = param1, param2

        # Input filenames
        rN = 'c_{0}_a_{1}_e_{2}_truth.txt'.format(c,a,e)
        gN = 'c_{0}_a_{1}_e_{2}_genotypes.txt'.format(c,a,e)
        rN, gN = [os.path.join(folder,x) for x in [rN, gN]]
        
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

        # | actual alpha - predicted alpha |
        alpha_error = abs(float(true_alpha[0]) - float(predict_alpha[0]))

        true_var = [k for k in true_cancer_dict.keys() if true_cancer_dict[k]!=true_normal_dict[k]]
        true_non_var = [k for k in true_cancer_dict.keys() if true_cancer_dict[k]==true_normal_dict[k]]
    
        assert len(true_var) + len(true_non_var) == len(true_cancer_dict.keys())

        predict_var = [k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]!=predict_normal_dict[k]]
        predict_non_var = [k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]==predict_normal_dict[k]]
        
        assert len(predict_var) + len(predict_non_var) == len(predict_cancer_dict.keys())

        assert tp + fp == len(predict_var)
        assert tn + fn == len(predict_non_var)

#        assert tp + fn == len(predict_var)
#        assert fp + tn == len(true_non_var)
        
        # try:
        #     assert fp + tn == len(predict_non_var)
        # except:
        #     print 'fp', fp
        #     print 'tp', tp
        #     print 'fn', fn
        #     print 'tn', tn
        #     print 'true var', len(true_var)
        #     print 'true non var', len(true_non_var)
        #     print 'predict var', len(predict_var)
        #     print 'predict non var', len(predict_non_var)
        #     0 / asdf

        z = False
        if z:
            print
            print a,c,e
            print rN
            print gN
            print folder, alpha_error, true_alpha, predict_alpha
            print 'fp', fp
            print 'tp', tp
            print 'fn', fn
            print 'tn', tn
            print 'true var', len(true_var)
            print 'true non var', len(true_non_var)
            print 'predict var', len(predict_var)
            print 'predict non var', len(predict_non_var)
            print true_var                                       

        tp2 = tp / float(tp + fp) if (tp+fp != 0) else None
        fp2 = fp / float(tp + fp) if (tp+fp != 0) else None
        tn2 = tn / float(tn + fn) if (tn+fn != 0) else None
        fn2 = fn / float(tn + fn) if (tn+fn != 0) else None

#        tp,fp,tn,fn=tp2,fp2,tn2,fn2

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

    tp_file.flush()
    fp_file.flush()
    fn_file.flush()
    tn_file.flush()
    alpha_file.flush()

tp_file.close()
fp_file.close()
fn_file.close()
tn_file.close()
alpha_file.close()

tv_file.close()
tnv_file.close()
pv_file.close()
pnv_file.close()
