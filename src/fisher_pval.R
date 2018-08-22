library(metap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

cell_marker = args[1]
tail = args[2]
input_folder = args[3]
siglim = as.numeric(args[4])
### extract filenames of the cell marker
file_list = list.files(input_folder, pattern=paste('^', cell_marker, '(.*)', tail, '$', sep='') )
print(file_list)
### read files of the cell marker
data_matrix = NULL
for (file in file_list){
	d = read.table(paste(input_folder, file, sep=''), header = F)
	d[d>308] = 308
	d[d<0] = 0
	data_matrix = cbind(data_matrix, d[,])
}
### get fisher method combined p-value
get_fisher_p = function(x){
	#print(x)
	if (length(x)!=1){
		x_p = 10^(-x)
		fp = sumlog(x_p)$p
		if (fp<=0.1^siglim){
			fp = 0.1^siglim
		}		
	} else{
		fp = 10^(-x)
		if (fp<=0.1^siglim){
			fp = 0.1^siglim
		}		
	}

	fp_neglog10 = -log10(fp)
	return(fp_neglog10)
}

fisher_p = apply(data_matrix, 1, function(x) get_fisher_p(x))

### write output
write.table(fisher_p, paste(cell_marker, '.fisher_p.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### write FRiP & SNRs
pk_sig = fisher_p[fisher_p>-log10(0.001)]
bg_sig = fisher_p[fisher_p<=-log10(0.001)]
FRiP = sum(pk_sig) / sum(bg_sig)
SNR = mean(pk_sig) / mean(bg_sig)

write.table(c(FRiP, SNR), paste(cell_marker, '.fisher_p.frip_snr.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

