#library(RColorBrewer)

method_list = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
ct_list = c('LSK_BM', 'HPC7', 'CMP', 'MEP', 'CFU_E_ad', 'G1E', 'ER4', 'ERY_ad', 'ERY_fl', 'CFUMK', 'MK_imm_ad', 'GMP', 'MONO_BM', 'NEU', 'NK_SPL', 'B_SPL', 'T_CD4_SPL', 'T_CD8_SPL')
t = 2
ct_list_plotting = c('LSK', 'HPC7', 'CMP', 'MEP', 'CFUE', 'G1E', 'ER4', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'GMP', 'MON', 'NEU', 'NK', 'B', 'TCD4', 'TCD8')




pkn_list_mat = c()
pkr_list_mat = c()

for (m in method_list){
pkn_list = c()
pkr_list = c()
for (c in ct_list){
### get input data
filename = paste(c, '_atac_macs_pk/', c, '.atac.meanrc.', m, '.macs2pk.ubg.NBQ.', t, '.txt', sep='')
print(filename)
d_tmp = read.table(filename, skip=1, header=F)
### get pkn & pkr
pkn_tmp = dim(d_tmp)[1]
pkr_tmp = sum(d_tmp[,3]-d_tmp[,2])
### append pkn & pkr
pkn_list = c(pkn_list, pkn_tmp)
pkr_list = c(pkr_list, pkr_tmp)
}
### get pkn & pkr matrix 
pkn_list_mat = cbind(pkn_list_mat, pkn_list)
pkr_list_mat = cbind(pkr_list_mat, pkr_list)
}

### get colnames & rownames
rownames(pkn_list_mat) = ct_list_plotting
colnames(pkn_list_mat) = method_list
rownames(pkr_list_mat) = ct_list_plotting
colnames(pkr_list_mat) = method_list

### get ylim
pkn_ylim = c(0, max(pkn_list_mat)+1000)
pkr_ylim = c(0, max(pkr_list_mat)+1000)

### get method colors
#method_colors = rev(brewer.pal(8,"Set3")[c(1,2,3,4,5)])
method_colors = c('gray', 'seagreen1', 'plum1', 'orange1', 'cadetblue1')
method_colors = c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1')

### plotting function
get_mergeplot = function(sig_mat, ct_names, outputname, method_colors, ylim_all){
pdf(outputname, width=6, height=6)
par(mar=c(5.5,4.1,2.1,2.1))
sig_vec = sig_mat[,1]
plot(1:length(sig_vec), sig_vec, pch = 16, col = method_colors[1], ylim=ylim_all, xaxt="n", xlab='', ylab='')
lines(1:length(sig_vec), sig_vec, lty=1, col=method_colors[1], lwd=2)
points(1:length(sig_vec), sig_vec, pch=1, col='black')
axis(1, at=1:length(sig_vec), labels=ct_names, col.axis="black", las=2, cex.axis=1.5)
for (i in 2:dim(sig_mat)[2]){
sig_vec = sig_mat[,i]
points(1:length(sig_vec), sig_vec, pch=16, col=method_colors[i])
lines(1:length(sig_vec), sig_vec, lty=1, col=method_colors[i], lwd=2)
points(1:length(sig_vec), sig_vec, pch=1, col='black')
}
dev.off()
}

### plotting
get_mergeplot(pkn_list_mat, ct_list_plotting, paste('barplot_pkn_pkr/atac.meanrc.allmethods.macs2pk.ubg.NBQ.', t, '.pkn.pdf', sep=''), method_colors, pkn_ylim)
get_mergeplot(pkr_list_mat, ct_list_plotting, paste('barplot_pkn_pkr/atac.meanrc.allmethods.macs2pk.ubg.NBQ.', t, '.pkr.pdf', sep=''), method_colors, pkr_ylim)



#get_barplot(t(pkn_list_mat), ct_list_plotting, paste('barplot_pkn_pkr/atac.meanrc.allmethods.macs2pk.NBQ.', t, '.pkn.pdf', sep=''), method_colors, pkn_ylim)
### plot barplot
#for (i in 1:dim(method_list)[2]){
### plot barplot 
#get_barplot(pkn_list_mat[,i], ct_list_plotting, paste('barplot_pkn_pkr/atac.meanrc.', method_list[i], '.macs2pk.NBQ.', t, '.pkn.pdf', sep=''), pkn_ylim)
#get_barplot(pkr_list_mat[,i], ct_list_plotting, paste('barplot_pkn_pkr/atac.meanrc.', method_list[i], '.macs2pk.NBQ.', t, '.pkr.pdf', sep=''), pkr_ylim)
#}




