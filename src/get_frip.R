### get parameters
args = commandArgs(trailingOnly=TRUE)
input_mat = args[1]
output = args[2]

### read input
all_sig = read.table(input_mat, header=F)

### get z
get_z = function(x){
	xz = (x - mean(x)) / sd(x)
	return(xz)
}

### get fdr
get_fdr = function(x){
	z = get_z(x)
	zp = pnorm(-abs(z))
	zpfdr = p.adjust(zp)
	return(zpfdr)
}

### get frip
get_frip = function(x, zpfdr, thresh) {
	pk_binary = (zpfdr < thresh)
	pks = sum(x[pk_binary])
	all = sum(x)
	frip = pks / all
	return(frip)
}

frip_all = c()
for (i in 1:dim(all_sig)[2]){
	print(i)
	sig_tmp = all_sig[,i]
	print('get fdr')
	sig_tmp_fdr = get_fdr(sig_tmp)
	print('get frip')
	sig_tmp_frip = get_frip(sig_tmp, sig_tmp_fdr, 0.05)
	frip_all = c(frip_all, sig_tmp_frip)
}

write.table(frip_all, output, sep='\t', quote=F, col.names=F, row.names=F)





