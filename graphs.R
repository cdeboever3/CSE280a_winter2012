library(gplots)

parse.table <- function(filename) {
  return(as.matrix(read.table(filename, check.names=FALSE)))
}

tp = parse.table('Figures/coverage_vs_alpha.tp.txt')
fp = parse.table('Figures/coverage_vs_alpha.fp.txt')
fn = parse.table('Figures/coverage_vs_alpha.fn.txt')
tn = parse.table('Figures/coverage_vs_alpha.tn.txt')

tp.rate = tp / (tp + fn)
fp.rate = fp / (tn + fp)
  
alpha = parse.table('Figures/coverage_vs_alpha.alpha.txt')

files <- list(list(table=tp.rate, filename='Figures/coverage_vs_alpha.tp_rate', title='True Positive Rate'),
              list(table=fp.rate, filename='Figures/coverage_vs_alpha.fp_rate', title='False Positive Rate'),
              list(table=fp, filename='Figures/coverage_vs_alpha.fp', title='# of False Positives'),
              list(table=fn, filename='Figures/coverage_vs_alpha.fn', title='# of False Negatives'),
              list(table=alpha, filename='Figures/coverage_vs_alpha.alpha', title='Alpha Error\n = |Actual alpha - Predicted Alpha|'))


for (i in seq(along=files)) {
  graph.table = files[[i]][['table']]
  graph.filename = files[[i]][['filename']]
  graph.title = files[[i]][['title']]
  
  png(paste(graph.filename,'.png',sep=''))
  heatmap.2(graph.table, Rowv=NA, Colv=NA, dendrogram="none", xlab='alpha (fraction of normal genome)', ylab='Average Coverage', hline=c(), vline=c(), trace="none")
  title(graph.title)
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
