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
mk_all = apply(file_list, 1, function(x) unlist(strsplit(x[1], '[.]'))[2])
mk_uniq = unique(mk_all)

for (mk in mk_uniq){
	print(mk)
	mk_binary = c()
	file_list_mk = file_list[mk_all==mk,1]
	###
	for (file in file_list_mk){
		print(file)
		signal_tmp = scan(file)
		signal_tmp_fdr = get_fdr(signal_tmp)
		signal_tmp_binary = (signal_tmp_fdr < thresh) * 1
		mk_binary_all = cbind(mk_binary, signal_tmp_binary)
	}
	mk_binary_all_pk = (apply(mk_binary_all, 1, prod) > 0 ) * 1
	mk_binary_all_bg = (apply(mk_binary_all, 1, sum) == 0 ) * 1
	output_mk_cpk = paste(output, '.', mk, '.cpk.txt', sep='')
	output_mk_cbg = paste(output, '.', mk, '.cbg.txt', sep='')
	write.table(mk_binary_all_pk, output_mk_cpk, sep='\t', quote=F, col.names=F, row.names=F)
	write.table(mk_binary_all_bg, output_mk_cbg, sep='\t', quote=F, col.names=F, row.names=F)
}


