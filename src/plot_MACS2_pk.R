library(RColorBrewer)

### methods & cell type list
ct_list = c('LSK_BM', 'HPC7', 'CMP', 'MEP', 'CFU_E_ad', 'G1E', 'ER4', 'ERY_fl', 'ERY_ad', 'CFUMK', 'MK_imm_ad', 'GMP', 'MONO_BM', 'NEU', 'NK_SPL', 'B_SPL', 'T_CD4_SPL', 'T_CD8_SPL')
methods_vec = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
methods_vec = c('raw', 'TSnorm', 'QTnorm', 'S3norm')

pkn_mat = c()
pkl_mat = c()
for (method in methods_vec){
### get pkn and pkl
pkn_vec = c()
pkl_vec = c()
for (ct in ct_list){
filename = paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.', method, '.macs2pk.ppois.10.txt', sep='')
print(filename)
pk_tmp = read.table(filename, skip=1, header=F)
pkn = dim(pk_tmp)[1]
pkl = sum(pk_tmp[,3]-pk_tmp[,2])
pkn_vec = c(pkn_vec, pkn)
pkl_vec = c(pkl_vec, pkl)
}
pkn_mat = rbind(pkn_mat, pkn_vec)
pkl_mat = rbind(pkl_mat, pkl_vec)
}

colnames(pkn_mat) = ct_list
rownames(pkn_mat) = methods_vec
colnames(pkl_mat) = ct_list
rownames(pkl_mat) = methods_vec


colors = rev(brewer.pal(8,"Set3")[1:5])[c(1,2,4,5)]

ylim_ud = c(0, max(pkn_mat))
pdf('pkn_pkl/pkn.pdf', width=4, height=4)
plot(1:dim(pkn_mat)[2], pkn_mat[1,], pch=16, col=colors[1], ylim = ylim_ud, xaxt="n", xlab='', ylab='')
lines(1:dim(pkn_mat)[2], pkn_mat[1,], col=colors[1])
for (i in 2:dim(pkn_mat)[1]){
points(1:dim(pkn_mat)[2], pkn_mat[i,], pch=16, col=colors[i])
lines(1:dim(pkn_mat)[2], pkn_mat[i,], col=colors[i])
}
axis(1, at=1:dim(pkn_mat)[2],labels=ct_list, col.axis="black", las=2, cex.axis=1)
dev.off()



ylim_ud = c(min(pkl_mat), max(pkl_mat))
pdf('pkn_pkl/pkl.pdf')
plot(1:dim(pkl_mat)[2], pkl_mat[1,], pch=16, col=colors[1], ylim = ylim_ud)
lines(1:dim(pkl_mat)[2], pkl_mat[1,], col=colors[1])
for (i in 2:dim(pkl_mat)[1]){
points(1:dim(pkl_mat)[2], pkl_mat[i,], pch=16, col=colors[i])
lines(1:dim(pkl_mat)[2], pkl_mat[i,], col=colors[i])
}
dev.off()


