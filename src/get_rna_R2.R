library(edgeR)
library(LSD)
library(ggplot2)
library(e1071) 

### read RNA seq
d0 = read.table('rnaHtseqCountsall.merged.pc.txt', header=FALSE)
d0_sig = d0[,-c(1)]
### add colnames
ct_list = c('CFU_E_ad', 'CMP', 'ERY_ad', 'GMP', 'MK_imm_ad', 'LSK_BM', 'MEP', 'MONO_BM', 'NEU', 'ER4', 'G1E')
colnames(d0_sig) = ct_list

### read gene coordinates
regions = read.table('gencode.vM4.annotation.pc.sorted.bed', header=FALSE)
regions_length = regions[,3]-regions[,2]

### edgeR normalize rna-seq
dgList <- DGEList(counts=d0_sig, genes=d0[,1])
dgList <- calcNormFactors(dgList, method="TMM")

### get RPKM
dgList_rpkm0 = rpkm(dgList, gene.length = regions_length)
#dgList_rpkm0 = cpm(dgList)
png('test.png')
heatscatter(dgList_rpkm0[,1],dgList_rpkm0[,6], log='xy')
abline(0,1,col='red')
dev.off()

### get rpkm threshold to filtering noise
dgList_rpkm_max = apply(dgList_rpkm0, 1, max)
rpkm_lim=-1.5
pdf('hist_dgList_rpkm_max.pdf', width=4, height=4)
hist(log2(dgList_rpkm_max), breaks=50)
abline(v=rpkm_lim, col='red', lwd=1.5, lty=2)
box()
dev.off()

### filtering genes
used_id_dgList_rpkm = (log2(dgList_rpkm_max+0.1)>rpkm_lim) >0
print(sum(used_id_dgList_rpkm))

### scale RNA-seq by all cell type mean and sd
dgList_rpkm = dgList_rpkm0
dgList_rpkm = as.matrix(log2(dgList_rpkm+0.1))
#dgList_rpkm = (dgList_rpkm-mean(dgList_rpkm))/sd(dgList_rpkm)


###############################
### get signal ct dif function
get_sigdif = function(x, i, j){
	c1sig = x[,i]
	c2sig = x[,j]
	cdif = c1sig-c2sig
	return(cdif)
}
### get R2 function
getR2 = function(obs, pred){
	R2_output = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
	return(R2_output)
}
### matrix to 2c table
get_frip_2c = function(x, name_end){
	x_2c = data.frame()
	for (i in 1:dim(x)[2]){
		x_2c = rbind(x_2c, cbind(x[,i], rep(paste(colnames(x)[i], name_end, sep=''), dim(x)[1])))
	}
	x_2c[,1] = as.numeric(as.character(x_2c[,1]))
	colnames(x_2c) = c('R2', 'Method')
	return(x_2c)
}
### matrix to 2c table
plot_heatscatter = function(obs, pred, outputname, lim, logxy, mainname){
	pdf(outputname, width=5, height=5)
	heatscatter(pred, obs, xlim=lim, ylim=lim, log=logxy, main=mainname)
	abline(0,1, lwd=1.5)
	abline(v=0, lwd=1.5)
	abline(h=0, lwd=1.5)
	dev.off()
}

###############################
methods_list = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')

R2o_mat = c()
R2i_mat = c()

set.seed(2018)
train_id = sample(dim(dgList_rpkm)[1], round(dim(dgList_rpkm)[1]*0.9, 0))

for (m in methods_list){
print(m)
### read h3k4me3 signal
h3k4me3 = read.table(paste('h3k4me3.TSSexp5kb.', m, '.pc.txt', sep=''), header=F)
h3k4me3_sig0 = h3k4me3[,-1]
h3k4me3_sig = as.matrix(log2(h3k4me3_sig0+0.1))
### output performace
R2o = c()
R2i = c()
for (i in 1:11){
print(i)
### get ct_i data
rna_i_train = (dgList_rpkm[train_id,i])
rna_i_test = (dgList_rpkm[-train_id,i])
h3k4me3_i_train = h3k4me3_sig[train_id,i]
h3k4me3_i_test = h3k4me3_sig[-train_id,i]
### svr training & testing
#lm_fit = lm(rna_i_train~h3k4me3_i_train)
#rna_i_pred = h3k4me3_i_test * lm_fit$coefficients[2] + lm_fit$coefficients[1]
rna_h3k4me3_train_i = as.data.frame(cbind(rna_i_train, h3k4me3_i_train))
colnames(rna_h3k4me3_train_i) = c('rna', 'h3k4me3')
svr_fit = svm(rna~h3k4me3, data=rna_h3k4me3_train_i)
h3k4me3_i_test = as.data.frame(h3k4me3_i_test)
colnames(h3k4me3_i_test) = c('h3k4me3')
rna_i_pred = predict(svr_fit, newdata=h3k4me3_i_test)
### get performance
ct_R2i = getR2(rna_i_test, rna_i_pred)
### append performance vector
R2i = c(R2i, ct_R2i)
### plot heatscatterplot
plot_output_names = paste('heatscatter/', ct_list[i], '_', ct_list[i], '.', m, '.hscatter.pdf', sep='')
plot_heatscatter(rna_i_test, rna_i_pred, plot_output_names, c(-5, 8), '', R2i)
for (j in 1:11){
if (i<j){
print(j)
### get ct j data
rna_j_test = (dgList_rpkm[-train_id,j])
h3k4me3_j_test = h3k4me3_sig[-train_id,j]
### predict using trained model
#rna_j_pred = h3k4me3_j_test * lm_fit$coefficients[2] + lm_fit$coefficients[1]
h3k4me3_j_test = as.data.frame(h3k4me3_j_test)
colnames(h3k4me3_j_test) = c('h3k4me3')
rna_j_pred = predict(svr_fit, newdata=h3k4me3_j_test)
### get performance
ct_R2o = getR2(rna_j_test, rna_j_pred)
### append performance vector
R2o = c(R2o, ct_R2o)
### plot heatscatterplot
plot_output_names = paste('heatscatter/', ct_list[i], '_', ct_list[j], '.', m, '.hscatter.pdf', sep='')
plot_heatscatter(rna_j_test, rna_j_pred, plot_output_names, c(-5, 8), '', R2o)
}
}
}
R2o_mat = cbind(R2o_mat, R2o)
R2i_mat = cbind(R2i_mat, R2i)
}

### add column names
colnames(R2o_mat) = methods_list
colnames(R2i_mat) = methods_list

R2o_mat_2c = get_frip_2c(R2o_mat, '_o')
R2i_mat_2c = get_frip_2c(R2i_mat, '_i')
R2_mat_2c = rbind(R2i_mat_2c, R2o_mat_2c)

pdf(paste('rna_vs_h3k4me3.pdf', sep=''), width=5.5, height=4)
p = ggplot(data = R2_mat_2c, aes(x=Method, y=(R2)))
p = p + geom_boxplot(aes(fill = Method))
p = p + geom_point(aes(y=(R2), group=Method), position = position_dodge(width=0.1))
p = p + scale_fill_manual(values=rep(c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1', 'gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1'),each = 1))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_x_discrete(breaks=unique(R2_mat_2c[,2]), labels=unique(R2_mat_2c[,2]))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12))
p = p + geom_hline(yintercept = 0, colour="gray", linetype="dashed", lwd=0.8)
p = p + geom_vline(xintercept = 5.5, colour="gray", linetype="dashed", lwd=0.8)
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
p = p + ylim(-1, 1)
plot(p)
dev.off()


png('test.png')
heatscatter(cdif_rna, cdif_h3k4me3)
abline(0,1)
abline(v=0)
abline(h=0)
dev.off()














methods_check = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
methods_check_name = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
methods_check_shuf = c('raw_shuf', 'TSnorm_shuf', 'MAnorm_shuf', 'QTnorm_shuf', 'S3norm_shuf')
methods_mix = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm', 'raw_shuf', 'TSnorm_shuf', 'MAnorm_shuf', 'QTnorm_shuf', 'S3norm_shuf')

r_mat = c()
r_mat_sp = c()
r_mat_shuf = c()
p_mat = c()
p_mat_shuf = c()
p_mat_ratio = c()
p_mat_shuf_ratio = c()
ct_pair_vec = c()


set.seed(2018)
shuffle_id1 = sample(1:11,11,replace=F)
shuffle_id2 = sample(1:11,11,replace=F)

for (i in c(9)){
	for (j in c(1:11)){
		if (i<j){
			k=0
			ct_pair_vec = rbind(ct_pair_vec, paste(ct_list[i], '-', ct_list[j], sep=''))
			r_vec = c()
			r_vec_sp = c()
			r_vec_shuf = c()
			p_vec = c()
			p_vec_shuf = c()
			p_vec_ratio = c()
			p_vec_shuf_ratio = c()

			m = 'rcnorm'
			d_raw0 = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+0.1)
			d_raw0 = (d_raw0 - mean(d_raw0))/sd(d_raw0)
			rna_dif0 = dgList_rpkm[used_id_dgList_rpkm,i]-dgList_rpkm[used_id_dgList_rpkm,j]
			hist_dif0 = d_raw0[used_id_dgList_rpkm,i]-d_raw0[used_id_dgList_rpkm,j]

			rss0 = (sum((rna_dif0 - hist_dif0)^2))
			for (m in methods_check){
				k = k+1
				png(paste('tss_h3k4me3.pcsorted.check.', m, '.', ct_list[i], '_', ct_list[j], '.png', sep=''))
				if (m != 'rcznorm'){
					d_raw = log10(as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))+0.1)
				} else {
					d_raw = as.matrix(read.table(paste('tss_h3k4me3.pcsorted.', m, '.txt', sep=''), header=F))
				}
				d_raw = (d_raw - mean(d_raw))/sd(d_raw)
				rna_dif = dgList_rpkm[used_id_dgList_rpkm,i]-dgList_rpkm[used_id_dgList_rpkm,j]
				hist_dif = d_raw[used_id_dgList_rpkm,i]-d_raw[used_id_dgList_rpkm,j]
				hist_dif_shuffle = d_raw[used_id_dgList_rpkm,shuffle_id1[i]]-d_raw[used_id_dgList_rpkm,shuffle_id1[j]]
				perform = 1 - sum((rna_dif - hist_dif)^2) / sum((rna_dif - mean(rna_dif))^2)
				#perform = sum((rna_dif>=0) == (hist_dif>=0))
				#perform = mean((rna_dif - hist_dif)^2)
				perform_shuf = 1 - sum((rna_dif - hist_dif_shuffle)^2) / sum((rna_dif - mean(rna_dif))^2)
				perform_ratio = (sum((rna_dif - hist_dif)^2)) / rss0
				perform_shuf_ratio = (sum((rna_dif - hist_dif_shuffle)^2)) / rss0
				cor0 = cor(rna_dif, hist_dif, method = 'pearson')
				cor0_sp = cor(rna_dif, hist_dif, method = 'spearman')
				cor0_shuf = cor(rna_dif, hist_dif_shuffle, method = 'pearson')
				r_vec[k] = cor0
				r_vec_sp[k] = cor0_sp
				r_vec_shuf[k] = cor0_shuf
				p_vec[k] = perform
				p_vec_shuf[k] = perform_shuf
				p_vec_ratio[k] = perform_ratio
				p_vec_shuf_ratio[k] = perform_shuf_ratio
				#heatscatter(dgList_rpkm[used_id_dgList_rpkm,i]-dgList_rpkm[used_id_dgList_rpkm,j], d_raw[used_id_dgList_rpkm,i]-d_raw[used_id_dgList_rpkm,j], log='', main=paste('R = ', toString(round(cor0, digits=3)), '; Ratio of RSS = ', toString(round(perform, digits=3)), sep=''), xlim=c(-3,3), ylim=c(-3,3) )
				plot(dgList_rpkm[used_id_dgList_rpkm,i]-dgList_rpkm[used_id_dgList_rpkm,j], d_raw[used_id_dgList_rpkm,i]-d_raw[used_id_dgList_rpkm,j], pch=16, log='', main=paste('R = ', toString(round(cor0, digits=3)), '; Ratio of RSS = ', toString(round(perform, digits=3)), sep=''), xlim=c(-5,5), ylim=c(-5,5) )
				abline(0,1,col='red')
				abline(v=0,col='blue')
				abline(h=0,col='blue')
				print(paste(ct_list[i], '_', ct_list[j], ': ', m, ': R = ', toString(round(cor0, digits=3))) )
				print(paste(ct_list[i], '_', ct_list[j], ': ', m, ': RSS = ', toString(round(perform, digits=3))) )
				dev.off()
			}
			r_mat = rbind(r_mat, r_vec)
			r_mat_sp = rbind(r_mat_sp, r_vec_sp)
			r_mat_shuf = rbind(r_mat_shuf, r_vec_shuf)
			p_mat = rbind(p_mat, p_vec)
			p_mat_shuf = rbind(p_mat_shuf, p_vec_shuf)
			p_mat_ratio = rbind(p_mat_ratio, p_vec_ratio)
			p_mat_shuf_ratio = rbind(p_mat_shuf_ratio, p_vec_shuf_ratio)
		}
	}
}



colnames(r_mat) = methods_check
colnames(r_mat_sp) = methods_check
colnames(r_mat_shuf) = methods_check
colnames(p_mat) = methods_check
colnames(p_mat_shuf) = methods_check
colnames(p_mat_ratio) = methods_check
colnames(p_mat_shuf_ratio) = methods_check

rownames(r_mat) = ct_pair_vec
rownames(r_mat_sp) = ct_pair_vec
rownames(r_mat_shuf) = ct_pair_vec
rownames(p_mat) = ct_pair_vec
rownames(p_mat_shuf) = ct_pair_vec
rownames(p_mat_ratio) = ct_pair_vec
rownames(p_mat_shuf_ratio) = ct_pair_vec

write.table(r_mat, 'r_mat.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(r_mat_sp, 'r_mat_sp.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )

write.table(r_mat_shuf, 'r_mat_shuf.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat, 'p_mat.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat_shuf, 'p_mat_shuf.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat_ratio, 'p_mat_ratio.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )
write.table(p_mat_shuf_ratio, 'p_mat_shuf_ratio.txt', quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t' )

r_mat = read.table('r_mat.txt', header = TRUE)
r_mat_sp = read.table('r_mat_sp.txt', header = TRUE)

r_mat_shuf = read.table('r_mat_shuf.txt', header = TRUE)
p_mat = read.table('p_mat.txt', header = TRUE)
p_mat_shuf = read.table('p_mat_shuf.txt', header = TRUE)
p_mat_ratio = read.table('p_mat_ratio.txt', header = TRUE)
p_mat_shuf_ratio = read.table('p_mat_shuf_ratio.txt', header = TRUE)


library(seewave)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

cor_mat_for_boxplot = c()
cor_mat_for_boxplot_sp = c()
cor_mat_for_boxplot_shuf = c()
p_mat_for_boxplot = c()
p_mat_for_boxplot_shuf = c()
p_mat_for_boxplot_ratio = c()
p_mat_for_boxplot_shuf_ratio = c()

for (i in c(1:length(methods_check))){
	### r_mat
	r_mat_tmp = r_mat[,i]
	method_tmp = rep(methods_check_name[i], length(r_mat[,i]))
	cor_mat_for_boxplot = rbind(cor_mat_for_boxplot, cbind(method_tmp, r_mat_tmp))
	### r_mat
	r_mat_tmp = r_mat_sp[,i]
	method_tmp = rep(methods_check_name[i], length(r_mat_sp[,i]))
	cor_mat_for_boxplot_sp = rbind(cor_mat_for_boxplot_sp, cbind(method_tmp, r_mat_tmp))
	### r_mat_shuf
	r_mat_tmp = r_mat_shuf[,i]
	method_tmp = rep(paste(methods_check_name[i], '_shuf', sep=''), length(r_mat[,i]))
	cor_mat_for_boxplot_shuf = rbind(cor_mat_for_boxplot_shuf, cbind(method_tmp, r_mat_tmp))
	### p_mat
	p_mat_tmp = p_mat[,i]
	method_tmp = rep(methods_check_name[i], length(p_mat[,i]))
	p_mat_for_boxplot = rbind(p_mat_for_boxplot, cbind(method_tmp, p_mat_tmp))
	### p_mat_shuf
	p_mat_tmp_shuf = p_mat_shuf[,i]
	method_tmp = rep(paste(methods_check_name[i], '_shuf', sep=''), length(p_mat[,i]))
	p_mat_for_boxplot_shuf = rbind(p_mat_for_boxplot_shuf, cbind(method_tmp, p_mat_tmp_shuf))
	### p_mat
	p_mat_tmp_ratio = p_mat_ratio[,i]
	method_tmp_ratio = rep(methods_check_name[i], length(p_mat_ratio[,i]))
	p_mat_for_boxplot_ratio = rbind(p_mat_for_boxplot_ratio, cbind(method_tmp, p_mat_tmp_ratio))
	### p_mat_shuf
	p_mat_tmp_shuf_ratio = p_mat_shuf_ratio[,i]
	method_tmp = rep(paste(methods_check_name[i], '_shuf', sep=''), length(p_mat_ratio[,i]))
	p_mat_for_boxplot_shuf_ratio = rbind(p_mat_for_boxplot_shuf_ratio, cbind(method_tmp, p_mat_tmp_shuf_ratio))
}

cor_mat_for_boxplot = as.data.frame(cor_mat_for_boxplot)
colnames(cor_mat_for_boxplot) = c('method', 'R')
cor_mat_for_boxplot$method = factor(cor_mat_for_boxplot$method, levels = methods_check_name,ordered = TRUE)
cor_mat_for_boxplot[,2] = apply(cor_mat_for_boxplot, 1, function(x) as.numeric(x[2]))
pdf(paste('R', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_name))) 
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

cor_mat_for_boxplot_sp = as.data.frame(cor_mat_for_boxplot_sp)
colnames(cor_mat_for_boxplot_sp) = c('method', 'R')
cor_mat_for_boxplot_sp$method = factor(cor_mat_for_boxplot_sp$method, levels = methods_check_name,ordered = TRUE)
cor_mat_for_boxplot_sp[,2] = apply(cor_mat_for_boxplot_sp, 1, function(x) as.numeric(x[2]))
pdf(paste('R_sp', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot_sp, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_name))) 
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

cor_mat_for_boxplot_shuf = as.data.frame(cor_mat_for_boxplot_shuf)
colnames(cor_mat_for_boxplot_shuf) = c('method', 'R')
cor_mat_for_boxplot_shuf$method = factor(cor_mat_for_boxplot_shuf$method, levels = methods_check_shuf,ordered = TRUE)
cor_mat_for_boxplot_shuf[,2] = apply(cor_mat_for_boxplot_shuf, 1, function(x) as.numeric(x[2]))
pdf(paste('R.shuf', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = cor_mat_for_boxplot_shuf, aes(x=method, y=R)) 
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=R, group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check)))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

p_mat_for_boxplot_od = p_mat_for_boxplot
p_mat_for_boxplot = as.data.frame(p_mat_for_boxplot)
colnames(p_mat_for_boxplot) = c('method', 'Performance')
p_mat_for_boxplot$method = factor(p_mat_for_boxplot$method, levels = methods_check_name,ordered = TRUE)
p_mat_for_boxplot[,2] = apply(p_mat_for_boxplot, 1, function(x) as.numeric(x[2]))
pdf(paste('performance', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot, aes(x=method, y=(Performance)) )
#p = p + ylim(c(-5,1))
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
#p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

p_mat_for_boxplot_ratio = as.data.frame(p_mat_for_boxplot_ratio)
colnames(p_mat_for_boxplot_ratio) = c('method', 'Performance')
p_mat_for_boxplot_ratio$method = factor(p_mat_for_boxplot_ratio$method, levels = methods_check_name,ordered = TRUE)
p_mat_for_boxplot_ratio[,2] = apply(p_mat_for_boxplot_ratio, 1, function(x) as.numeric(x[2]))
pdf(paste('performance_ratio', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_ratio, aes(x=method, y=(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


p_mat_for_boxplot_shuf = as.data.frame(p_mat_for_boxplot_shuf)
colnames(p_mat_for_boxplot_shuf) = c('method', 'Performance')
p_mat_for_boxplot_shuf$method = factor(p_mat_for_boxplot_shuf$method, levels = methods_check_shuf,ordered = TRUE)
p_mat_for_boxplot_shuf[,2] = apply(p_mat_for_boxplot_shuf, 1, function(x) as.numeric(x[2]))
pdf(paste('performance.shuf', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_shuf, aes(x=method, y=(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_shuf))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()

p_mat_for_boxplot_shuf_ratio = as.data.frame(p_mat_for_boxplot_shuf_ratio)
colnames(p_mat_for_boxplot_shuf_ratio) = c('method', 'Performance')
p_mat_for_boxplot_shuf_ratio$method = factor(p_mat_for_boxplot_shuf_ratio$method, levels = methods_check_shuf,ordered = TRUE)
p_mat_for_boxplot_shuf_ratio[,2] = apply(p_mat_for_boxplot_shuf_ratio, 1, function(x) as.numeric(x[2]))
pdf(paste('performance.shuf', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_shuf_ratio, aes(x=method, y=(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
#p = p + scale_fill_manual(values=rep("deepskyblue", length(methods_check_shuf))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


p_mat_for_boxplot_all = as.data.frame(rbind(p_mat_for_boxplot, p_mat_for_boxplot_shuf))
colnames(p_mat_for_boxplot_all) = c('method', 'Performance')
p_mat_for_boxplot_all$method = factor(p_mat_for_boxplot_all$method, levels = methods_mix,ordered = TRUE)
p_mat_for_boxplot_all[,2] = apply(p_mat_for_boxplot_all, 1, function(x) as.numeric(x[2]))
pdf(paste('performance.all', '.box.pdf', sep=''))#, width=500, height=500)
p = ggplot(data = p_mat_for_boxplot_all, aes(x=method, y=1/(Performance)) )
p = p + geom_boxplot(aes(fill = method))
p = p + geom_point(aes(y=1/(Performance), group=method), position = position_dodge(width=0.75))
p = p + geom_hline(yintercept = 1, color="black", linetype="dashed")
p = p + scale_fill_manual(values=rep(c("deepskyblue", 'gray'), length(methods_mix))) + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()


