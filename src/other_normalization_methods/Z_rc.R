args = commandArgs(trailingOnly=TRUE)
ref = args[1]
target = args[2]
target_norm = args[3]

small_num = 0
target_sig  = scan(target)
#ref_sig = scan(ref)
target_sig_Z = (target_sig - mean(target_sig))/sd(target_sig)
write.table(target_sig_Z, target_norm, quote=F, sep='\t', col.names=F, row.names=F)

