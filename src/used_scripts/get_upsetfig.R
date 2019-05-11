library(UpSetR)

### get parameters
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]

#input = 'G1E_atac_macs_pk/G1E.atac.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt'
#output = 'upset/G1E.atac.meanrc.macs2pk.NBQ.2.pooled.pk.pdf'

data = read.table(input, header = FALSE)

method_list = c('raw', 'TSnorm', 'MAnorm', 'QTnorm', 'S3norm')
method_colors = c('gray', 'seagreen1', 'plum1', 'orange1', 'dodgerblue1')

#set.seed(2018)
#used_id = sample(dim(data)[1], 10000)
#data_s = as.data.frame(data[used_id,-1])

data_binary = data[,-c(1,2,3)]
colnames(data_binary) = method_list

data_binary[data_binary>0] = 1

pdf(output, height=3, width=7)
upset(data_binary, nsets = 10, keep.order = TRUE, nintersects = NA, empty.intersections=0, group.by = 'degree', sets = method_list,
 scale.intersections="identity", scale.sets="identity", order.by = c('degree'),decreasing = c(TRUE, TRUE),
 sets.bar.color = method_colors)
dev.off()


