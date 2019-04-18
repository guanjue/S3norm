
method_vec = c('raw', 'TSnorm', 'S3norm', 'QTnorm')#, 'POIS', 'NBP', 'S3norm_NBP')
color_vec = c('black', 'blue', 'red', 'orange', 'green', 'purple', 'cyan')
thresh_vec = seq(0,100,2)

pk_info_all = c()
pk_info = c()
for (i in thresh_vec){
	print(i)
	filename = paste('CMP.atac.meanrc.', method_vec[1], '.macs2pk.qpois.', i,'.txt', sep='')
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
	filename = paste('CMP.atac.meanrc.', m, '.macs2pk.qpois.', i,'.txt', sep='')
	d = read.table(filename, skip=1)
	pk_len = sum(d[,3]-d[,2])
	#pk_num = dim(d)[1]
	pk_info = c(pk_info, pk_len)
}
pk_info_all = cbind(pk_info_all, pk_info)
}


win = 3
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


pdf('CMP.atac.meanrc.macs2pk.qpois.thresh_curve.pdf', height=21)
par(mfrow=c(3,1))
plot(pk_info_all[,1], pk_info_all[,2], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all[,-1]), max(pk_info_all[,-1])))
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all[,2+k], col = color_vec[k+1])
}
### plot 1d
plot(pk_info_all[,1], pk_info_all_1d_all[,1], log='y', type="l", col = color_vec[1], ylim=c(min(pk_info_all_1d_all), max(pk_info_all_1d_all)))
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all_1d_all[,1+k], col = color_vec[k+1])
}
### plot 2d
plot(pk_info_all[,1], pk_info_all_2d_all[,1], log='', type="l", col = color_vec[1], ylim=c(min(pk_info_all_2d_all), max(pk_info_all_2d_all)))
for (k in 1:length(method_vec[-1])){
lines(pk_info_all[,1], pk_info_all_2d_all[,1+k], col = color_vec[k+1])
}
dev.off()


pdf('CMP.atac.meanrc.macs2pk.qpois.1dhist.pdf', height=30)
par(mfrow=c(dim(pk_info_all_1d_all)[2]+1,1))
for (i in 1:dim(pk_info_all_1d_all)[2]){
len_thresh = dim(pk_info_all_1d_all)[1]
hist(pk_info_all_1d_all[-c(1:win, (len_thresh-2):len_thresh),i], breaks=50)
}
hist(pk_info_all_1d_all, breaks=50)
dev.off()









