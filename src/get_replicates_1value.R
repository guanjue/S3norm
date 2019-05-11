library(ggplot2)
library(LSD)

### matrix to 2c table
get_frip_2c = function(x, name_end){
	x_2c = data.frame()
	for (i in 1:dim(x)[2]){
		x_2c = rbind(x_2c, cbind(x[,i], rep(paste(colnames(x)[i], name_end, sep=''), dim(x)[1])))
	}
	x_2c[,1] = as.numeric(as.character(x_2c[,1]))
	colnames(x_2c) = c('R', 'Method')
	return(x_2c)
}

### plot 
plot_heatscatter = function(x, outputname){
	png(outputname)
	lim = c(min(x)-0.1, max(x)+0.1)
	heatscatter(x[,1], x[,2], xlim=lim, ylim=lim)
	abline(0, 1, lwd=1.5)
	dev.off()
}


### add colnames
ct_list = c('CFU_E_ad', 'CFUMK', 'CMP', 'ERY_ad', 'GMP', 'MK_imm_ad', 'LSK_BM', 'MEP', 'NEU', 'ER4', 'G1E')
#ct_list = c('CFU_E_ad', 'CFUMK', 'CMP', 'ERY_ad', 'ER4', 'G1E')

methods_list = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
#methods_list = c('raw', 'QTnorm', 'S3norm')

R_mat_2c = c()

for (m in methods_list){
print(m)
R_mat = c()
for (ct in ct_list){
print(ct)
filename = paste(ct, 'rep12.atac.meanrc.', m, '.bw.tab.txt', sep='')
d = read.table(filename, header=F)
plot_heatscatter(log2(d+1), paste('replicate_cor/', filename, '.png', sep=''))
dm = rowMeans(d)
R_vec = cor(log2(d[,1]+1), log2(d[,2]+1))
R_mat = c(R_mat, R_vec)
dat = cbind(log2(d[,1]+1), log2(d[,2]+1))
#R_vec = est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
#pk_idr = R_vec$IDR
#R_mat = c(R_mat, mean(pk_idr))
}
output_name = paste('R_mat.', m, '.R.pr.txt', sep='')
write.table(R_mat, output_name, quote=F, col.names=F, row.names=F)
###
R_mat_2c = cbind(R_mat_2c, R_mat)
}


rownames(R_mat_2c) = ct_list
colnames(R_mat_2c) = c(methods_list)
R_mat_2c = get_frip_2c(R_mat_2c, '')
#R_mat_2c[,1] = -log10(R_mat_2c[,1])

pdf(paste('replicate_cor/replicates_R.pdf', sep=''), width=5.5, height=4)
p = ggplot(data = R_mat_2c, aes(x=Method, y=(R)))
p = p + geom_boxplot(aes(fill = Method))
p = p + geom_point(aes(y=(R), group=Method), position = position_dodge(width=0.1))
p = p + scale_fill_manual(values=rep(c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1'), each=1))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_x_discrete(breaks=unique(R_mat_2c[,2]), labels=unique(R_mat_2c[,2]))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12))
#p = p + geom_hline(yintercept = 0, colour="gray", linetype="dashed", lwd=0.8)
#p = p + geom_vline(xintercept = 5.5, colour="gray", linetype="dashed", lwd=0.8)
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
p = p + ylim(0, 1)
p = p + theme(legend.position = "none")
plot(p)
dev.off()


set.seed(2018)
used_id = sample(dim(d)[1], 10000)
png('replicate_cor/test.png')
plot(d[used_id,1], d[used_id,2], log='xy', xlim=c(1,10000), ylim=c(1,10000))
abline(0,1)
dev.off()



