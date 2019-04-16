args = commandArgs(trailingOnly=TRUE)
ref = args[1]
target = args[2]
target_norm = args[3]

target_sig  = scan(target)
ref_sig = scan(ref)

target_sig_qtnorm = ref_sig[order(ref_sig)][rank(target_sig)]

write.table(target_sig_qtnorm, target_norm, quote=F, sep='\t', col.names=F, row.names=F)

