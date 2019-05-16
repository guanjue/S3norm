colSD = function(x){
	return(apply(x, 2, sd))
}

R2 = function(x1, x2){
	shuf_id = sample(length(x2), length(x2))
	x2_shuf = x2[shuf_id]
	R2i = 1-sum((x1-x2)^2) / sum((x1-mean(x1))^2)
	return(R2i)
}

mse = function(x1, x2){
	return(mean((x1-x2)^2))
}

ct1 = 'G1E'
method_vec = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
method_color_vec = c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1')

setwd(paste('/gpfs/scratch/gzx103/vision/all_final_data/rc_norm/sigdif_pk/', ct1, '_GMP_pooledctpk', sep=''))

#method_vec = c('raw', 'QTnorm', 'S3norm')
db = read.table(paste('../../', ct1, '_atac_macs_pk/', ct1, '.atac.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.binary.txt', sep=''), header=F)
db_b = db[,-c(1:3)]
db_b[db_b>1] = 1
db_b_l = apply(db_b, 1, function(x) paste(x, collapse=''))

ct_vec = c('CFU_E_ad', 'CFUMK', 'CMP', 'ER4', 'ERY_ad', 'GMP', 'LSK_BM', 'MEP', 'NEU')
m_ID = 0
cor_rep_in_mat = c()
for (m in method_vec){
m_ID = m_ID+1
print(m)
### read data
d0 = read.table(paste('', ct1, '.GMP.repsig.', m, '.pooledct.txt', sep=''), header=F)
d = d0
d[,2:5] = log2(d[,2:5]+1)
dm = rowMeans(d[,c(2,4)])
dm_range = max(dm)-min(dm)
t1 = 1
n = 10
dm_range_d = dm_range/n
interval = 1/n
t_vec = 1-rev(seq(0, 1, interval))
### get cor vec
cor_rep_in = c()
cor_rep_shuf_mat = c()
#print(t_vec)
#print(t_vec+interval)
###
set.seed(2018)
t_id = 0
for (t2 in t_vec){
t_id = t_id+1
d_rep1_tmp = (d[(dm<=quantile(dm,1)) & (dm>quantile(dm,t2)),2])
d_rep2_tmp = (d[(dm<=quantile(dm,1)) & (dm>quantile(dm,t2)),4])
#d_rep1_tmp = (d[(dm<=max(dm)-(t_id-1)*dm_range_d) & (dm>max(dm)-(t_id-0)*dm_range_d),2])
#d_rep2_tmp = (d[(dm<=max(dm)-(t_id-1)*dm_range_d) & (dm>max(dm)-(t_id-0)*dm_range_d),4])
#cor_rep_in = c(cor_rep_in, cor(d_rep1_tmp, d_rep2_tmp))
cor_rep_in = c(cor_rep_in, cor(d_rep1_tmp,d_rep2_tmp))
set.seed(2018)
cor_rep_shuf_vec = c()
for (si in 1:1){
used_id = sample(length(d_rep2_tmp), length(d_rep2_tmp))
#cor_rep_shuf_vec = c(cor_rep_shuf_vec, mean(as.matrix((d_rep1_tmp-d_rep2_tmp[used_id])^2)))
cor_rep_shuf_vec = c(cor_rep_shuf_vec, cor(d_rep1_tmp,d_rep2_tmp[used_id]))	
}
cor_rep_shuf_mat = cbind(cor_rep_shuf_mat, cor_rep_shuf_vec)
}

cor_rep_in_mat = rbind(cor_rep_in_mat, cor_rep_in)


###plot
pdf(paste('test.', m, '.pdf', sep=''))
plot(100*(1-t_vec), cor_rep_in, ylim=c(-2,1), type='l', col=method_color_vec[m_ID])
lines(100*(1-t_vec), colMeans(cor_rep_shuf_mat))
points(100*(1-t_vec), colMeans(cor_rep_shuf_mat), pch=16, col='gray30')
lines(100*(1-t_vec), colMeans(cor_rep_shuf_mat)+1.96*colSD(cor_rep_shuf_mat), lty=2, col='gray30')
lines(100*(1-t_vec), colMeans(cor_rep_shuf_mat)-1.96*colSD(cor_rep_shuf_mat), lty=2, col='gray30')
points(100*(1-t_vec), cor_rep_in, pch=16, col=method_color_vec[m_ID])
abline(v=sum(db_b[,1]==1)/dim(db_b)[1]*100, col=method_color_vec[1])
abline(v=sum(db_b[,2]==1)/dim(db_b)[1]*100, col=method_color_vec[2])
abline(v=sum(db_b[,3]==1)/dim(db_b)[1]*100, col=method_color_vec[3])
abline(v=sum(db_b[,4]==1)/dim(db_b)[1]*100, col=method_color_vec[4])
abline(v=sum(db_b[,5]==1)/dim(db_b)[1]*100, col=method_color_vec[5])
dev.off()
###
png(paste('test.', m, '.png', sep=''))
plot(d[,2], d[,4], log='', pch=16, xlim=log2(c(1,2000)), ylim=log2(c(1,2000)))
points(d[dm>quantile(dm,0.9),2], d[dm>quantile(dm,0.9),4], col='red')
points(d[dm>quantile(dm,0.95),2], d[dm>quantile(dm,0.95),4], col='blue')
points(d[dm>quantile(dm,0.975),2], d[dm>quantile(dm,0.975),4], col='green')
points(d[dm>quantile(dm,0.99),2], d[dm>quantile(dm,0.99),4], col='orange')
abline(0,1)
dev.off()
png(paste('test.o.', m, '.png', sep=''))
plot(d[,2], d[,3], log='', pch=16)
points(d[dm>quantile(dm,0.9),2], d[dm>quantile(dm,0.9),3], col='red')
points(d[dm>quantile(dm,0.95),2], d[dm>quantile(dm,0.95),3], col='blue')
points(d[dm>quantile(dm,0.975),2], d[dm>quantile(dm,0.975),3], col='green')
points(d[dm>quantile(dm,0.99),2], d[dm>quantile(dm,0.99),3], col='orange')
abline(0,1)
dev.off()
}

###plot
pdf(paste('test.', 'allmethods', '.pdf', sep=''))
plot(100*(1-t_vec), cor_rep_in_mat[1,], ylim=c(0,1), type='l', col=method_color_vec[1], log='')
points(100*(1-t_vec), cor_rep_in_mat[1,], pch=16, col=method_color_vec[1])
for (i in 2:5){
lines(100*(1-t_vec), cor_rep_in_mat[i,], col=method_color_vec[i])
points(100*(1-t_vec), cor_rep_in_mat[i,], pch=16, col=method_color_vec[i])
}
abline(v=sum(db_b[,1]==1)/dim(db_b)[1]*100, col=method_color_vec[1])
abline(v=sum(db_b[,2]==1)/dim(db_b)[1]*100, col=method_color_vec[2])
abline(v=sum(db_b[,3]==1)/dim(db_b)[1]*100, col=method_color_vec[3])
abline(v=sum(db_b[,4]==1)/dim(db_b)[1]*100, col=method_color_vec[4])
abline(v=sum(db_b[,5]==1)/dim(db_b)[1]*100, col=method_color_vec[5])
dev.off()







