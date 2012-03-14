library(gplots)

fp = read.table('coverage_vs_alpha.fp.txt')
fp = as.matrix(fp)
png('False_Positive.png')
heatmap.2(fp, Rowv=NA, Colv=NA, labCol=(1:19)/20, xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace='none')
#title('False Positive Rate \n= (# of predicted somatic mutations that are actually germline) \n -------------/ (# of predicted somatic mutations)')
title('False Positive Rate')
dev.off()

print(fp)

fn = read.table('coverage_vs_alpha.fn.txt')
fn = as.matrix(fn)
png('False_Negative.png')
heatmap.2(fn, Rowv=NA, Colv=NA, labCol=(1:19)/20, xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace='none')
#title('False Negative Rate \n= (# of germline mutations that are actually somatic) / (# of germline mutations)')
title('False Negative Rate')
dev.off()

print(fn)

alpha = read.table('coverage_vs_alpha.alpha.txt')
alpha = as.matrix(alpha)
png('Alpha_Error.png')
heatmap.2(alpha, Rowv=NA, Colv=NA, labCol=(1:19)/20, xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace='none')
title('Alpha Error\n = |Actual alpha - Predicted alpha|')
dev.off()

print(alpha)
