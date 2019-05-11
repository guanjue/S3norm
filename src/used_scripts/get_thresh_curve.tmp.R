### get parameters
args = commandArgs(trailingOnly=TRUE)
ct = args[1]


#method_vec = c('raw', 'TSnorm', 'QTnorm', 'S3norm', 'POIS', 'NBP', 'S3norm_NBP')
method_vec = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')#, 'POIS', 'NBP', 'S3norm_NBP')
method_vec = c('raw', 'TSnorm', 'QTnorm', 'S3norm')
color_vec = c('black', 'blue', 'red', 'orange', 'green', 'purple', 'cyan')
thresh_vec = seq(0,100,2)
thresh_vec = seq(0,100,1)
thresh_vec = seq(0,40,1)

pk_info_all = c()
pk_info = c()
for (i in thresh_vec){
	print(i)
	filename = paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.', method_vec[1], '.macs2pk.qpois.', i,'.txt', sep='')
	d = read.table(filename, skip=1)
	pk_len = sum(d[,3]-d[,2])
	#pk_num = dim(d)[1]
	pk_info = rbind(pk_info, c(i, pk_len))
}
pk_info_all = cbind(pk_info_all, pk_info)
### other methods

for (m in method_vec[-1]){
pk_info = c()
for (i in thresh_vec){
	print(i)
	print(m)
	filename = paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.', m, '.macs2pk.qpois.', i,'.txt', sep='')
	d = read.table(filename, skip=1)
	pk_len = sum(d[,3]-d[,2])
	#pk_num = dim(d)[1]
	pk_info = c(pk_info, pk_len)
}
pk_info_all = cbind(pk_info_all, pk_info)
}


win = 2
pk_info_all_1d_all = c()
for (j in 2:dim(pk_info_all)[2]){
pk_info_all_1d = c()
for (i in (1+win):(dim(pk_info_all)[1]-win)){
	d = (log2(pk_info_all[i-win,j]) - log2(pk_info_all[i+win,j])) / (2*win)
	pk_info_all_1d = c(pk_info_all_1d, d)
}
pk_info_all_1d = c(rep(pk_info_all_1d[1], win), pk_info_all_1d, rep(pk_info_all_1d[length(pk_info_all_1d)], win))
pk_info_all_1d_all = cbind(pk_info_all_1d_all, pk_info_all_1d)
}


#win = 3
pk_info_all_2d_all = c()
for (j in 1:dim(pk_info_all_1d_all)[2]){
pk_info_all_2d = c()
for (i in (1+win):(dim(pk_info_all_1d_all)[1]-win)){
	d = log2(pk_info_all_1d_all[i-win,j]) - log2(pk_info_all_1d_all[i+win,j])
	pk_info_all_2d = c(pk_info_all_2d, d)
}
pk_info_all_2d = c(rep(pk_info_all_2d[1], win), pk_info_all_2d, rep(pk_info_all_2d[length(pk_info_all_2d)], win))
pk_info_all_2d_all = cbind(pk_info_all_2d_all, pk_info_all_2d)
}



#win = 3
pk_info_all_3d_all = c()
for (j in 1:dim(pk_info_all_2d_all)[2]){
pk_info_all_3d = c()
for (i in (1+win):(dim(pk_info_all_2d_all)[1]-win)){
	d = log2(pk_info_all_2d_all[i-win,j]) - log2(pk_info_all_2d_all[i+win,j])
	pk_info_all_3d = c(pk_info_all_3d, d)
}
pk_info_all_3d = c(rep(pk_info_all_3d[1], win), pk_info_all_3d, rep(pk_info_all_3d[length(pk_info_all_3d)], win))
pk_info_all_3d_all = cbind(pk_info_all_3d_all, pk_info_all_3d)
}



get_gamma_thresh = function(x, fdr_thresh){
derivative_vector = as.vector(x)
derivative_vector_rmtop = derivative_vector[derivative_vector<quantile(derivative_vector, 0.95)]
s = var(derivative_vector_rmtop) / mean(derivative_vector_rmtop)
a = mean(derivative_vector_rmtop)^2 / var(derivative_vector_rmtop)
gamma_p = (1-pgamma(as.vector(x), shape=a, scale=s))
#gamma_p_fdr = p.adjust(gamma_p, 'fdr')
gamma_limit = min(derivative_vector[gamma_p<fdr_thresh])
for (i in c(1:5)){
derivative_vector = derivative_vector[gamma_p>=fdr_thresh]
derivative_vector_rmtop = derivative_vector
s = var(derivative_vector_rmtop) / mean(derivative_vector_rmtop)
a = mean(derivative_vector_rmtop)^2 / var(derivative_vector_rmtop)
gamma_p = (1-pgamma(as.vector(x), shape=a, scale=s))
#gamma_p_fdr = p.adjust(gamma_p, 'fdr')
gamma_limit = min(derivative_vector[gamma_p<fdr_thresh])
if (sum(gamma_p<fdr_thresh)==0){
break
}
}
return(gamma_limit)
}



gamma_thresh_vec = c()
pdf(paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.macs2pk.qpois.1dhist.pdf', sep=''), height=20)
par(mfrow=c(dim(pk_info_all_1d_all)[2],1))
for (i in 1:dim(pk_info_all_1d_all)[2]){
len_thresh = dim(pk_info_all_1d_all)[1]
derivative_1d = pk_info_all_1d_all[-c(1:win, (len_thresh-2):len_thresh),i]
gamma_thresh = get_gamma_thresh(derivative_1d, 0.05)
gamma_thresh_vec = c(gamma_thresh_vec, gamma_thresh)
hist(derivative_1d, breaks=100)
abline(v=gamma_thresh)
box()
}
dev.off()



pdf(paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.macs2pk.qpois.thresh_curve.pdf', sep=''), height=14)
par(mfrow=c(2,1))
### plot 1d
plot(pk_info_all[,1], pk_info_all_1d_all[,1], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all_1d_all), max(pk_info_all_1d_all)))
derivative_lim = max(pk_info_all_1d_all[,1][(pk_info_all_1d_all[,1]<gamma_thresh_vec[1])])
derivative_lim_thresh = pk_info_all[,1][which(pk_info_all_1d_all[,1]==derivative_lim)]
points(derivative_lim_thresh, derivative_lim, col=color_vec[1])
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all_1d_all[,1+k], col = color_vec[k+1])
derivative_lim = max(pk_info_all_1d_all[,k+1][(pk_info_all_1d_all[,k+1]<gamma_thresh_vec[k+1])])
derivative_lim_thresh = pk_info_all[,1][which(pk_info_all_1d_all[,k+1]==derivative_lim)]
points(derivative_lim_thresh, derivative_lim, col=color_vec[k+1])
}
### plot pk num
region_lim_vec = c()
region_lim_thresh_vec = c()
plot(pk_info_all[,1], pk_info_all[,2], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all[,-1]), max(pk_info_all[,-1])))
region_lim = max(pk_info_all[,2][(pk_info_all_1d_all[,1]<gamma_thresh_vec[1])])
region_lim_thresh = pk_info_all[,1][which(pk_info_all[,2]==region_lim)]
points(region_lim_thresh, region_lim, col=color_vec[1])
###
region_lim_vec = c(region_lim_vec, region_lim)
region_lim_thresh_vec = c(region_lim_thresh_vec, region_lim_thresh)
###
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all[,2+k], col = color_vec[k+1])
region_lim = max(pk_info_all[,2+k][(pk_info_all_1d_all[,k+1]<gamma_thresh_vec[k+1])])
region_lim_thresh = pk_info_all[,1][which(pk_info_all[,2+k]==region_lim)]
points(region_lim_thresh, region_lim, col=color_vec[k+1])
### 
region_lim_vec = c(region_lim_vec, region_lim)
region_lim_thresh_vec = c(region_lim_thresh_vec, region_lim_thresh)
}
write.table(cbind(region_lim_thresh_vec, region_lim_vec), paste(ct, '_atac_macs_pk/', ct, '.atac.lim.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
dev.off()




gamma_thresh_vec = c()
pdf(paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.macs2pk.qpois.counthist.pdf', sep=''), height=20)
par(mfrow=c(dim(pk_info_all)[2]-1,1))
for (i in 2:dim(pk_info_all)[2]){
len_thresh = dim(pk_info_all)[1]
pk_info_all_log2 = log2(pk_info_all[,i])
gamma_thresh = get_gamma_thresh(pk_info_all_log2, 0.05)
gamma_thresh_vec = c(gamma_thresh_vec, gamma_thresh)
hist(pk_info_all_log2, breaks=30)
abline(v=gamma_thresh)
box()
}
dev.off()


pdf(paste(ct, '_atac_macs_pk/', ct, '.atac.meanrc.macs2pk.qpois.c.thresh_curve.pdf', sep=''), height=14)
par(mfrow=c(2,1))
### plot 1d
plot(pk_info_all[,1], pk_info_all_1d_all[,1], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all_1d_all), max(pk_info_all_1d_all)))
derivative_lim = max(pk_info_all_1d_all[,1][(log2(pk_info_all[,2])<gamma_thresh_vec[1])])
derivative_lim_thresh = pk_info_all[,1][which(pk_info_all_1d_all[,1]==derivative_lim)]
points(derivative_lim_thresh, derivative_lim, col=color_vec[1])
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all_1d_all[,1+k], col = color_vec[k+1])
derivative_lim = max(pk_info_all_1d_all[,k+1][(log2(pk_info_all[,2+k])<gamma_thresh_vec[k+1])])
derivative_lim_thresh = pk_info_all[,1][which(pk_info_all_1d_all[,k+1]==derivative_lim)]
points(derivative_lim_thresh, derivative_lim, col=color_vec[k+1])
}
### plot pk num
region_lim_vec = c()
region_lim_thresh_vec = c()
pk_info_all = pk_info_all[-1,]
plot(pk_info_all[,1], pk_info_all[,2], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all[,-1]), max(pk_info_all[,-1])))
region_lim = max(pk_info_all[,2][(log2(pk_info_all[,2])<gamma_thresh_vec[1])])
region_lim_thresh = pk_info_all[,1][which(pk_info_all[,2]==region_lim)]
points(region_lim_thresh, region_lim, col=color_vec[1])
###
region_lim_vec = c(region_lim_vec, region_lim)
region_lim_thresh_vec = c(region_lim_thresh_vec, region_lim_thresh)
###
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all[,2+k], col = color_vec[k+1])
region_lim = max(pk_info_all[,2+k][(log2(pk_info_all[,2+k])<gamma_thresh_vec[k+1])])
region_lim_thresh = pk_info_all[,1][which(pk_info_all[,2+k]==region_lim)]
points(region_lim_thresh, region_lim, col=color_vec[k+1])
### 
region_lim_vec = c(region_lim_vec, region_lim)
region_lim_thresh_vec = c(region_lim_thresh_vec, region_lim_thresh)
}
write.table(cbind(region_lim_thresh_vec, region_lim_vec), paste(ct, '_atac_macs_pk/', ct, '.atac.lim.count.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
dev.off()



