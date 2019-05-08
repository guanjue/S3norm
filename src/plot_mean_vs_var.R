library(ggplot2)
library(RColorBrewer)

file_list = read.table('rc_list.atac.txt', header=F)
method_list = c('RAW', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')


get_2c = function(x, name_end){
        x_2c = data.frame()
        for (i in 1:dim(x)[2]){
                x_2c = rbind(x_2c, cbind(x[,i], rep(paste(colnames(x)[i], name_end, sep=''), dim(x)[1])))
        }
        x_2c[,1] = as.numeric(as.character(x_2c[,1]))
        colnames(x_2c) = c('Info', 'T')
        return(x_2c)
}

mean_var_mat_all = c()
for (j in 1:length(method_list)){
mean_vec = c()
var_vec = c()
ct_vec = c()
### get all info
for (i in 1:dim(file_list)[1]){
mean_file = read.table(paste(file_list[i,3], '.', file_list[i,4], '.pkbg.txt', sep=''))[j,2]
sd_file = read.table(paste(file_list[i,3], '.', file_list[i,4], '.pkbg.sd.txt', sep=''))[j,2]
mean_vec = c(mean_vec, mean_file)
var_vec = c(var_vec, sd_file^2)
ct_vec = c(ct_vec, as.character(file_list[i,3]))
}
### get mean var mat
mean_var_mat = cbind(mean_vec, var_vec)
rownames(mean_var_mat) = ct_vec
colnames(mean_var_mat) = c(paste(method_list[j], '_Mean', sep=''), paste(method_list[j], '_Var', sep=''))
mean_var_mat = (mean_var_mat)
mean_var_mat_all = cbind(mean_var_mat_all, mean_var_mat)
}

mean_var_mat_all = cbind(mean_var_mat_all[,c(1,3,5,7,9)], mean_var_mat_all[,c(2,4,6,8,10)])
mean_var_mat_all_2c = get_2c(mean_var_mat_all, '')

### plot boxplot
pdf(paste('Allmethods', '.mean_vs_var.boxplot.pdf', sep=''), width=5.5, height=4)
p = ggplot(data = mean_var_mat_all_2c, aes(x=T, y=log2(Info)))
p = p + geom_boxplot(aes(fill = T))
p = p + geom_point(aes(y=log2(Info), group=T), position = position_dodge(width=0.1))
p = p + scale_fill_manual(values=rep(rev(brewer.pal(8,"Set3")[1:5]), 2))
#p = p + scale_fill_brewer(palette="YlGnBu", direction=-1)
p = p + theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
p = p + scale_x_discrete(breaks=unique(mean_var_mat_all_2c[,2]), labels=unique(mean_var_mat_all_2c[,2]))
p = p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
p = p + ylim(0, 10.5)
plot(p)
dev.off()



