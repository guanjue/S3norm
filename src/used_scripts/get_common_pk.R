### get parameters
args = commandArgs(trailingOnly=TRUE)
file_list_file = args[1]
output = args[2]
thresh = as.numeric(args[3])

### get z
get_z = function(x){
	x_notop = x[x<=quantile(x, 0.99)]
        xz = (x - mean(x_notop)) / sd(x_notop)
        return(xz)
}

### get fdr
get_fdr = function(x){
        z = get_z(x)
        zp = pnorm(-abs(z))
        zpfdr = p.adjust(zp)
        return(zpfdr)
}


### read input
file_list = read.table(file_list_file, header=F)

common_pk = c()
for (i in 1:dim(file_list)[1]){
	print(file_list[i,3:4])
	sig_tmp = scan(paste(file_list[i,3], '.', file_list[i,4], '.meanrc.txt', sep=''))
	sig_tmp_fdr = get_fdr(sig_tmp)
	sig_tmp_fdr_pk = sig_tmp_fdr<thresh
	print(sum(sig_tmp_fdr_pk))
	common_pk = cbind(common_pk, sig_tmp_fdr_pk)
}


common_pk_binary = apply(common_pk, 1, prod)
common_bg_binary = apply(common_pk, 1, sum)  

write.table(common_pk_binary, paste(output,'.cpk.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
write.table(common_bg_binary, paste(output,'.cbg.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
