library(gplots)

parse.table <- function(filename) {
  return(as.matrix(read.table(filename, check.names=FALSE)))
}

folder = 'Figures_original/'

tp = parse.table(paste(folder,'alpha_vs_error.tp.txt',sep=''))
fp = parse.table(paste(folder,'alpha_vs_error.fp.txt',sep=''))
fn = parse.table(paste(folder,'alpha_vs_error.fn.txt',sep=''))
tn = parse.table(paste(folder,'alpha_vs_error.tn.txt',sep=''))

## tp.rate = tp / (tp + fn)
## fp.rate = fp / (tn + fp)
  
alpha = abs(parse.table(paste(folder,'alpha_vs_error.alpha.txt',sep='')))

tp = t(tp)
fp = t(fp)
fn = t(fn)
alpha = t(alpha)

print(fp)
print(tp)
print(alpha)

fn[1,1]=.001

files <- list(#list(table=tp.rate, filename=paste(folder,'alpha_vs_error_tp_rate',sep=''), title='True Positive Rate'),
              #list(table=fp.rate, filename=paste(folder,'alpha_vs_error_fp_rate',sep=''), title='False Positive Rate'),
              list(table=fp, filename=paste(folder,'alpha_vs_error_fp',sep=''), title='# of False Positives'),
              list(table=fn, filename=paste(folder,'alpha_vs_error_fn',sep=''), title='# of False Negatives'),
              list(table=alpha, filename=paste(folder,'alpha_vs_error_alpha',sep=''), title='Alpha Error\n = |Actual alpha - Predicted Alpha|'))

#alpha = abs(parse.table(paste(folder,'alpha_vs_error_improved.alpha.txt',sep='')))
#files <- list(list(table=alpha, filename=paste(folder,'alpha_vs_error_improved.alpha',sep=''),title='Alpha Error\n= |Actual alpha - Predicted alpha|'))

for (i in seq(along=files)) {
  graph.table = files[[i]][['table']]
  graph.filename = files[[i]][['filename']]
  graph.title = files[[i]][['title']]

  print(graph.table)
  png(paste(graph.filename,'.png',sep=''), width=480,height=480)
#  pdf(paste(graph.filename,'.png',sep=''))

  ## heatmap.2(graph.table, Rowv=NA, Colv=NA, dendrogram="none", xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace="none",key=TRUE, lmat=rbind( c(2,3,0,0),c(0,4,0,0),c(0,1,1,1)), lhei=c(0.5,1,3.0), lwid=c(0.5,2,2,2),cexRow=1.5,cexCol=1.5)

  heatmap.2(graph.table, Rowv=NA, Colv=NA, dendrogram="none", xlab='Seq Error Rate', ylab='alpha (fraction of normal genome)', hline=c(), vline=c(), trace="none",key=TRUE, lmat=rbind( c(2,3,0,0),c(0,4,0,0),c(0,1,1,1)), lhei=c(0.5,1,3.0), lwid=c(0.5,2,2,2),cexRow=1.5,cexCol=1.5)

  
#  title(main=graph.title, font.main=30)
  dev.off()
}  

## fp.name = 'Figures/coverage_vs_alpha.fp'
## fp = read.table(paste(fp.name,'.txt',sep=''), check.names=FALSE)
## fp = as.matrix(fp)
## png(paste(fp.name,'.png',sep=''))
## heatmap.2(fp, Rowv=NA, Colv=NA, dendrogram="none", xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace="none")
## title('False Positive Rate')
## dev.off()

## fn.name = 'Figures/coverage_vs_alpha.fn'
## fn = read.table(paste(fn.name,'.txt',sep=''), check.names=FALSE)
## fn = as.matrix(fn)
## png(paste(fn.name,'.png',sep=''))
## heatmap.2(fn, Rowv=NA, Colv=NA, dendrogram="none", xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace="none")
## title('False Negative Rate')
## dev.off()

## alpha.name = 'Figures/coverage_vs_alpha.alpha'
## alpha = read.table(paste(alpha.name,'.txt',sep=''), check.names=FALSE)
## alpha = as.matrix(alpha)
## png(paste(alpha.name,'.png',sep=''))
## heatmap.2(alpha, Rowv=NA, Colv=NA, dendrogram="none", xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace="none")
## title('Alpha Error\n = |Actual alpha - Predicted alpha|')
## dev.off()
