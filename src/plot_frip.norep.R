library(ggplot2)
library(RColorBrewer)

list = read.table('rc_list.atac.norep.txt', header=F)

ct_list = as.character(list[,3])
mk_list = as.character(list[,4])

methods = c('RAW', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')

merge_frips = c()

for (i in 1:length(ct_list)){
	ct_rep12merge = paste(ct_list[i], '.', mk_list[i], '.frip.txt', sep='')
	ct_rep12merge_a = as.matrix(read.table(ct_rep12merge, header=F))
	merge_frips = cbind(merge_frips, ct_rep12merge_a)
}


rownames(merge_frips) = methods
colnames(merge_frips) = ct_list
merge_frips = t(merge_frips)


get_frip_2c = function(x, name_end){
	x_2c = data.frame()
	for (i in 1:dim(x)[2]){
		x_2c = rbind(x_2c, cbind(x[,i], rep(paste(colnames(x)[i], name_end, sep=''), dim(x)[1])))
	}
	x_2c[,1] = as.numeric(as.character(x_2c[,1]))
	colnames(x_2c) = c('FRIP', 'Method')
	return(x_2c)
}

merge_frips_2c = get_frip_2c(merge_frips, '')

frip_2c = data.frame()
for (i in 1:length(methods)){
	used_id = (1+(i-1)*length(ct_list)):(length(ct_list)+(i-1)*length(ct_list))
	frip_2c = rbind(frip_2c, merge_frips_2c[used_id,])
}


frip_2c$Method = factor(frip_2c$Method, levels = unique(frip_2c$Method),ordered = TRUE)

pdf(paste('frip', '.norep.box.pdf', sep=''), width=8, height=4)
p = ggplot(data = frip_2c, aes(x=Method, y=FRIP)) 
p = p + geom_boxplot(aes(fill = Method))
p = p + geom_point(aes(y=FRIP, group=Method), position = position_dodge(width=0.75))
p = p + scale_fill_manual(values=rep(rev(brewer.pal(8,"YlGnBu")[1:5]),each = 1))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_x_discrete(breaks=unique(frip_2c[,2]), labels=unique(frip_2c[,2]))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
plot(p)
dev.off()



