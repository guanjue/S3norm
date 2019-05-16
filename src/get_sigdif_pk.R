args = commandArgs(trailingOnly=TRUE)
ct1 = args[1]
ct2 = args[2]

library(pheatmap)
library(mclust)
library(DESeq2)
library(LSD)

set.seed(2018)

### functions
colMedian = function(x){
xm = apply(x,2,median)
return(xm)
}


#ct1='CMP'
#ct2='G1E'
#m='S3norm'
#m='raw'
#m='QTnorm'
#m='MAnorm'

#CMP_G1E_sharedpk/CMP.G1E.repsig.S3norm.txt

#filename = paste(ct1, '_', ct2, '_sharedpk/', ct1, '.', ct2, '.repsig.', m, '.txt', sep='')

method_vec = c('raw', 'TSnorm', 'QTnorm', 'S3norm', 'MAnorm')
#method_vec = c('QTnorm', 'S3norm')

for (m in method_vec){

filename = paste('sigdif_pk/', ct1, '_', ct2, '_sharedpk/', ct1, '.', ct2, '.repsig.', m, '.txt', sep='')
d = read.table(filename, header=F)

### get colnames
colnames_name = c('ccRE_id', paste(ct1, '.rep1', sep=''), paste(ct2, '.rep1', sep=''), paste(ct1, '.rep2', sep=''), paste(ct2, '.rep2', sep=''))
colnames(d) = colnames_name
### get conditions
d_condition = cbind(c(ct1, ct2, ct1, ct2))
colnames(d_condition) = 'condition'

###### start DEseq2
### get DEseq object
dsig = d[,-1]
rownames(dsig) = d[,1]

### DEseq2
dds <- DESeqDataSetFromMatrix(countData = round(dsig,0), colData = d_condition, design= ~ condition)

### set reference
#dds$condition <- relevel(dds$condition, ref = "CMP")

### run DEseq2
dds <- DESeq(dds)

### set sizeFactors to 1
sizeFactors(dds) = rep(1.0, dim(dsig)[2])

### run DEseq2 without normalization
dds <- DESeq(dds)

### adjp thresh
adj_thresh = 0.05

### MAplot
res <- results(dds)

M = res$log2FoldChange
M[M>7]=7
M[M<(-7)]=-7
A = res$baseMean
png(paste('sigdif_pk/', ct1, '_', ct2, '_sharedpk/', ct1, '.', ct2, '.repsig.', m, '.png', sep=''), width=600, height=300)
par(mfrow=c(1,2))
plotMA(dds, alpha=adj_thresh, ylim=c(-7,7), xlim=c(1,1e+3))
heatscatter(A, M, log='x', ylim=c(-7,7), xlim=c(1,1e+3))
abline(h=0)
dev.off()

c1 = sum(((res$padj[!is.na(res$padj)])<adj_thresh) & (res$log2FoldChange[!is.na(res$padj)]>0))
c2 = sum(((res$padj[!is.na(res$padj)])<adj_thresh) & (res$log2FoldChange[!is.na(res$padj)]<0))

write.table(c(c1, c2), paste('sigdif_pk/', ct1, '_', ct2, '_sharedpk/', ct1, '.', ct2, '.repsig.', m, '.difpkcount.txt', sep=''))
print(m)
print(sum(((res$padj[!is.na(res$padj)])<adj_thresh) & (res$log2FoldChange[!is.na(res$padj)]>0)))
print(sum(((res$padj[!is.na(res$padj)])<adj_thresh) & (res$log2FoldChange[!is.na(res$padj)]<0)))

}


