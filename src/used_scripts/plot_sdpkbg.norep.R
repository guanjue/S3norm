library(ggplot2)
library(RColorBrewer)

list = read.table('rc_list.atac.txt', header=F)

ct_list = as.character(list[,3])
mk_list = as.character(list[,4])

methods = c('RAW', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')

merge_pk = c()
merge_bg = c()

for (i in 1:length(ct_list)){
	ct_rep12merge = paste(ct_list[i], '.', mk_list[i], '.pkbg.sd.txt', sep='')
	ct_rep12merge_a = as.matrix(read.table(ct_rep12merge, header=F))
	merge_pk = cbind(merge_pk, ct_rep12merge_a[,1])
	merge_bg = cbind(merge_bg, ct_rep12merge_a[,2])
}


rownames(merge_pk) = methods
colnames(merge_pk) = ct_list
merge_pk = t(merge_pk)
rownames(merge_bg) = methods
colnames(merge_bg) = ct_list
merge_bg = t(merge_bg)


get_frip_2c = function(x, name_end){
	x_2c = data.frame()
	for (i in 1:dim(x)[2]){
		x_2c = rbind(x_2c, cbind(x[,i], rep(paste(colnames(x)[i], name_end, sep=''), dim(x)[1])))
	}
	x_2c[,1] = as.numeric(as.character(x_2c[,1]))
	colnames(x_2c) = c('SD', 'Method_regions')
	return(x_2c)
}

merge_pk_2c = get_frip_2c(merge_pk, '_PK')
merge_bg_2c = get_frip_2c(merge_bg, '_BG')

frip_2c = data.frame()
for (i in 1:length(methods)){
	used_id = (1+(i-1)*length(ct_list)):(length(ct_list)+(i-1)*length(ct_list))
	frip_2c = rbind(frip_2c, merge_pk_2c[used_id,])
	frip_2c = rbind(frip_2c, merge_bg_2c[used_id,])
}

frip_2c$Method_regions = factor(frip_2c$Method_regions, levels = unique(frip_2c$Method_regions),ordered = TRUE)

pdf(paste('pkbg', '.atac.merged.box.sd.pdf', sep=''), width=5.5, height=4)
p = ggplot(data = frip_2c, aes(x=Method_regions, y=log2(SD))) 
p = p + geom_boxplot(aes(fill = Method_regions))
p = p + geom_point(aes(y=log2(SD), group=Method_regions), position = position_dodge(width=0.1))
p = p + scale_fill_manual(values=rep(c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1'),each = 2))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_x_discrete(breaks=unique(frip_2c[,2]), labels=unique(frip_2c[,2]))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
p = p + geom_hline(yintercept = log2(merge_pk[1,1]), colour="orange", linetype="dashed", lwd=0.8)
p = p + geom_hline(yintercept = log2(merge_bg[1,1]), colour="gray", linetype="dashed", lwd=0.8)
p = p + ylim(0, 9)
plot(p)
dev.off()




