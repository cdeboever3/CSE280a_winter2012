import os, glob

folder_list = ['coverage_vs_alpha', 'alpha_vs_error']

def get_param_1():
    pass

def get_param_2():
    pass

for folder in folder_list:
    
    print folder

    predict_files = glob.glob(os.path.join(folder, '*_genotypes.txt'))
    pileup_files = [x.split('_genotypes.txt')[0] + '.tsv' for x in predict_files]
    truth_files = [x.split('_genotypes.txt')[0] + '_truth.txt' for x in predict_files]


    # param_1 = 'coverage'
    # if folder=='coverage_vs_alpha':
    #     param_2 = 'alpha'
    # elif folder=='alpha_vs_error':
    #     param_2 = 'error'
    # print >>outF, '\t'.join([param_1, param_2, 'alphaError', 'tp', 'fp', 'tn', 'fn'])
    # print >>outF, '\t'.join([param_1, param_2, 'alphaError', 'tp', 'fp', 'tn', 'fn'])

    #for predict, pileup, truth in zip(predict_files, pileup_files, truth_files):

    covL = [2**i for i in range(2, 8)]
    alpha_step = 0.05
    alphaL = [round(alpha_step*i, 4) for i in range(1, int(1 / alpha_step))]

    errorL = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]

    if folder=='coverage_vs_alpha':
        param1_list = covL
        param2_list = alphaL
    else:
        param1_list = alphaL
        param2_list = errorL

    fp_file = open(folder + '.fp.txt', 'w')
    fn_file = open(folder + '.fn.txt', 'w')
    alpha_file = open(folder + '.alpha.txt', 'w')
    
    # Print headers
    print >>fp_file, '\t' + '\t'.join(map(str,param2_list))
    print >>fn_file, '\t' + '\t'.join(map(str,param2_list))
    print >>alpha_file, '\t' + '\t'.join(map(str,param2_list))
    
    for param1 in param1_list:

        fp_list, fn_list, alpha_list = [], [], []
        
        for param2 in param2_list:
            base = '%s_%s_%s_%s' % ('c' if folder=='coverage_vs_alpha' else 'a',
                                    param1,
                                    'a' if folder=='coverage_vs_alpha' else'e',
                                    param2)
            base = os.path.join(folder, base)
            truth = base + '_truth.txt'
            predict = base + '_genotypes.txt'
            pileup = base + '.tsv'

            # params = os.path.basename(truth).split('_truth.txt')[0].split('_')
            # param_1 = float(params[1 + params.index('a')])
            # if folder=='coverage_vs_alpha':
            #     param_2 = coverage= float(params[1 + params.index('c')])
            # elif folder=='alpha_vs_error':
            #     param_2 = error = float(params[1 + params.index('e')])

            # dict: genomic position --> actual genotype of cancer
            truth_lines = [x for x in open(truth).read().split('\n') if x!='']
            truth_pop_freq = [x.split('\t')[1:3] for x in truth_lines if x[0]=='#' and 'population frequencies' in x]
            assert len(truth_pop_freq)==1
            truth_pop_freq = truth_pop_freq[0]
            truth_lines = [x.split() for x in truth_lines if x[0]!='#']
            truth_cancer_dict = dict([(x[0]+'.'+x[1], set(x[3].upper())) for x in truth_lines])
            truth_normal_dict = dict([(x[0]+'.'+x[1], set(x[2].upper())) for x in truth_lines])

            # dict: genomic position --> predicted genotype of cancer
            predict_lines = [x for x in open(predict).read().split('\n') if x!='']
            predict_pop_freq = [x.split('\t')[1:3] for x in predict_lines if x[0]=='#' and 'population frequencies' in x]
            assert len(predict_pop_freq)==1        
            predict_pop_freq = predict_pop_freq[0]
            predict_lines = [x.split() for x in predict_lines if x[0]!='#']
            predict_cancer_dict = dict([(x[0]+'.'+x[1], set(x[3].upper())) for x in predict_lines])
            predict_normal_dict = dict([(x[0]+'.'+x[1], set(x[2].upper())) for x in predict_lines])
            a, b = truth_cancer_dict, truth_normal_dict
            c, d = predict_cancer_dict, predict_normal_dict

            fp = len([k for k in c.keys() if
                      c[k]!=d[k] and \
                          a.has_key(k) and \
                          a[k]==b[k] \
                      ])

            tp = len([k for k in c.keys() if
                      c[k]!=d[k] and \
                          a.has_key(k) and \
                          a[k]!=b[k] \
                          ])

            a, b = predict_cancer_dict, predict_normal_dict
            c, d = truth_cancer_dict, truth_normal_dict

            fn = len([k for k in c.keys() if
                      c[k]!=d[k] and \
                          a.has_key(k) and \
                          a[k]==b[k] \
                      ])

            tn = len([k for k in c.keys() if
                      c[k]==d[k] and \
                          a.has_key(k) and \
                          a[k]==b[k] \
                          ])

            alpha_error = abs(float(truth_pop_freq[0]) - float(predict_pop_freq[0]))

            truth_var = len([k for k in truth_cancer_dict.keys() if truth_cancer_dict[k]!=truth_normal_dict[k]])
            truth_non_var = len([k for k in truth_cancer_dict.keys() if truth_cancer_dict[k]==truth_normal_dict[k]])
            predict_var = len([k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]!=predict_normal_dict[k]])
            predict_non_var = len([k for k in predict_cancer_dict.keys() if predict_cancer_dict[k]==predict_normal_dict[k]])

            # print predict
            # print folder, coverage, alpha_error, truth_pop_freq, predict_pop_freq
            # print 'fp', fp
            # print 'tp', tp
            # print 'fn', fn
            # print 'tn', tn
            # print 'true var', truth_var
            # print 'true non var', truth_non_var
            # print 'predict var', predict_var
            # print 'predict non var', predict_non_var

            tp_rate = tp / float(tp + fp)
            fp_rate = fp / float(tp + fp)
            tn_rate = tn / float(tn + fn)
            fn_rate = fn / float(tn + fn)

            alpha_list.append(alpha_error)
            fp_list.append(fp_rate)
            fn_list.append(fn_rate)
#            fp_list.append(fp)
#            fn_list.append(fn)

        print >>fp_file, str(param1) + '\t' + '\t'.join(map(str,fp_list))
        print >>fn_file, str(param1) + '\t' + '\t'.join(map(str,fn_list))
        print >>alpha_file, str(param1) + '\t' + '\t'.join(map(str,alpha_list))

        fp_file.flush()
        fn_file.flush()
        alpha_file.flush()

    fp_file.close()
    fn_file.close()
    alpha_file.close()
