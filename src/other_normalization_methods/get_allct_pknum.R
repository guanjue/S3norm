### get parameters
args = commandArgs(trailingOnly=TRUE)
ct = args[1]
ct='CMP'
thresh = args[1]
thresh = 5

#method_vec = c('raw', 'TSnorm', 'QTnorm', 'S3norm', 'POIS', 'NBP', 'S3norm_NBP')
method_vec = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')#, 'POIS', 'NBP', 'S3norm_NBP')
#method_vec = c('raw', 'TSnorm')#, 'POIS', 'NBP', 'S3norm_NBP')
thresh_vec = seq(0,100,1)
thresh_vec = seq(2,30,1)
color_vec = c('black', 'blue', 'red', 'orange', 'green', 'purple', 'cyan')
ct_files = c("LSK_BM", "HPC7", "CMP", "MEP", "CFU_E_ad", "G1E", "ER4", "ERY_ad", "ERY_fl", "CFUMK", "MK_imm_ad", "GMP", "MONO_BM", "NEU", "NK_SPL", "B_SPL", "T_CD4_SPL", "T_CD8_SPL")


for (m in method_vec){
pk_info_all = c()
pk_info_pkn_all = c()
for (t in thresh_vec){
print(t)
pk_info_region = c()
pk_info_pk = c()
for (ct in ct_files){
print(ct)
filename = paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.', m, '.macs2pk.qpois.', t,'.txt', sep='')
	d = read.table(filename, skip=1)
	pk_len = sum(d[,3]-d[,2])
	pk_num = dim(d)[1]
	pk_info_region = rbind(pk_info_region, pk_len)
	pk_info_pk = rbind(pk_info_pk, pk_num)
}
pk_info_all = cbind(pk_info_all, pk_info_region)
pk_info_pkn_all = cbind(pk_info_pkn_all, pk_info_pk)
}

### get rownames & colnames
rownames(pk_info_all) = ct_files
rownames(pk_info_pkn_all) = ct_files
colnames(pk_info_all) = thresh_vec
colnames(pk_info_pkn_all) = thresh_vec

pk_info_all = t(pk_info_all)
pk_info_pkn_all = t(pk_info_pkn_all)

lim_used_region = c(min(pk_info_all), max(pk_info_all))
lim_used_pkn = c(min(pk_info_pkn_all), max(pk_info_pkn_all))

### pk region
#dir.create('pk_all_ct')
pdf(paste('pk_all_ct/pk_region.', m, '.pdf',sep=''), height=7)
par(mfrow=c(1,1))
ct_x = c(1:dim(pk_info_all)[2])
plot(ct_x, pk_info_all[1,], log='y', type="l", col = color_vec[1], ylim=lim_used_region)
points(ct_x, pk_info_all[1,], col=color_vec[1])
for (t in 1:(dim(pk_info_all)[1]-1)){
lines(ct_x, pk_info_all[t,], col = color_vec[1])
points(ct_x, pk_info_all[t,], col=color_vec[1])
}
dev.off()
### pkn
pdf(paste('pk_all_ct/pkn.', m, '.pdf',sep=''), height=7)
par(mfrow=c(1,1))
ct_x = c(1:dim(pk_info_pkn_all)[2])
plot(ct_x, pk_info_pkn_all[1,], log='y', type="l", col = color_vec[1], ylim=lim_used_pkn)
points(ct_x, pk_info_pkn_all[1,], col=color_vec[1])
for (t in 2:(dim(pk_info_pkn_all)[1])){
lines(ct_x, pk_info_pkn_all[t,], col = color_vec[1])
points(ct_x, pk_info_pkn_all[t,], col=color_vec[1])
}
dev.off()
}



for (t in thresh_vec){
pk_info_all_difm = c()
pk_info_pkn_all_difm = c()

for (m in method_vec){
pk_info_region = c()
pk_info_pk = c()
for (ct in ct_files){
print(ct)
filename = paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.', m, '.macs2pk.qpois.', t,'.txt', sep='')
	d = read.table(filename, skip=1)
	pk_len = sum(d[,3]-d[,2])
	pk_num = dim(d)[1]
	pk_info_region = rbind(pk_info_region, pk_len)
	pk_info_pk = rbind(pk_info_pk, pk_num)
}
pk_info_all_difm = cbind(pk_info_all_difm, pk_info_region)
pk_info_pkn_all_difm = cbind(pk_info_pkn_all_difm, pk_info_pk)
}


### get rownames & colnames
rownames(pk_info_all_difm) = ct_files
rownames(pk_info_pkn_all_difm) = ct_files
colnames(pk_info_all_difm) = method_vec
colnames(pk_info_pkn_all_difm) = method_vec

pk_info_all_difm = t(pk_info_all_difm)
pk_info_pkn_all_difm = t(pk_info_pkn_all_difm)

lim_used_region_difm = c(min(pk_info_pkn_all_difm), max(pk_info_pkn_all_difm))
lim_used_pkn_difm = c(min(pk_info_all_difm), max(pk_info_all_difm))

pdf(paste('pk_all_ct/pk_region.across.method.', t, '.pdf',sep=''), height=7)
par(mfrow=c(1,1))
ct_x = c(1:dim(pk_info_pkn_all_difm)[2])
plot(ct_x, pk_info_pkn_all_difm[1,], log='y', type="l", col = color_vec[1], ylim=lim_used_region_difm)
points(ct_x, pk_info_pkn_all_difm[1,], col=color_vec[1])
for (t in 2:(dim(pk_info_pkn_all_difm)[1])){
lines(ct_x, pk_info_pkn_all_difm[t,], col = color_vec[t])
points(ct_x, pk_info_pkn_all_difm[t,], col=color_vec[t])
}
dev.off()

pdf(paste('pk_all_ct/pkn.across.method.', t, '.pdf',sep=''), height=7)
par(mfrow=c(1,1))
ct_x = c(1:dim(pk_info_all_difm)[2])
plot(ct_x, pk_info_all_difm[1,], log='y', type="l", col = color_vec[1], ylim=lim_used_pkn_difm)
points(ct_x, pk_info_all_difm[1,], col=color_vec[1])
for (t in 2:(dim(pk_info_all_difm)[1])){
lines(ct_x, pk_info_all_difm[t,], col = color_vec[t])
points(ct_x, pk_info_all_difm[t,], col=color_vec[t])
}
dev.off()
}

