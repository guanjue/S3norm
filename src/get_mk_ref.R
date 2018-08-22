### get parameters
args = commandArgs(trailingOnly=TRUE)

parameters_file_list = args[1]
method = args[2]
output_name = args[3]
parameters_files = read.table(parameters_file_list, header=F)



FRiP_list = c()
SNR_list = c()
for (i in c(1:dim(parameters_files)[1])){
        print(i)
        print(parameters_files[i,1])
        file=toString(parameters_files[i,1])
        parameters = read.table(file, header=F)
        ### FRiP score is in the 5th column
        FRiP_list[i] = parameters[1,]
        SNR_list[i] = parameters[2,]
}
print(dim(parameters_files))
print(FRiP_list)
print(which.max(SNR_list))
print(SNR_list)
print(which.max(SNR_list))


if (method=='snr'){
	ref_file = sub('.frip_snr.txt', '.txt', toString(parameters_files[which.max(SNR_list),1]))
	frip_ref = FRiP_list[which.max(SNR_list)]
	SNR_ref = SNR_list[which.max(SNR_list)]
	write.table(c(ref_file, frip_ref, SNR_ref), output_name, sep='\t', quote=F, col.names=F, row.names=F)
	ref_file_final = sub('.fisher_p.frip_snr.txt', '.pknorm.ref.txt', toString(parameters_files[which.max(SNR_list),1]))
	parameters_files_name = apply(parameters_files, 1, function(x) sub('.frip_snr.txt', '.txt', x))
	pknorm_list = cbind(rep(ref_file_final, dim(parameters_files)[1]), parameters_files_name)
	write.table(pknorm_list, paste(output_name, '.info.txt', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
} else if (method=='frip'){
	ref_file = sub('.frip_snr.txt', '.txt', toString(parameters_files[which.max(FRiP_list),1]))
	frip_ref = FRiP_list[which.max(FRiP_list)]
	SNR_ref = SNR_list[which.max(FRiP_list)]
	write.table(c(ref_file, frip_ref, SNR_ref), output_name, sep='\t', quote=F, col.names=F, row.names=F)
	ref_file_final = sub('.fisher_p.frip_snr.txt', '.pknorm.ref.txt', toString(parameters_files[which.max(FRiP_list),1]))
	parameters_files_name = apply(parameters_files, 1, function(x) sub('.frip_snr.txt', '.txt', x))
	pknorm_list = cbind(rep(ref_file_final, dim(parameters_files)[1]), parameters_files_name)
	write.table(pknorm_list, paste(output_name, '.info.txt', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
}
