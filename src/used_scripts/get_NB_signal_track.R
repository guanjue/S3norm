### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
input_track_file = args[2]
cbg_track_file = args[3]
output_name = args[4]

#time Rscript get_NB_signal_track.R MK_imm_ad.atac.meanrc.S3norm.sort.txt MK_imm_ad.atac.meanrc.S3norm.sort.bg.bedgraph atac_commonpk.txt.cbg.sort.txt MK_imm_ad.atac.meanrc.S3norm.sort


### get 0-adjusted NB model
get_true_NB_prob_size = function(x){
	m=mean(x[x>0]);
	m2=mean(x[x>0]^2);
	p0 = length(which(x==0)) / length(x);
	p = m/(m2-m^2 * (1-p0));
	#if (p<0.1){
	#	p = 0.1
	#}
	#if (p>=0.9){
	#	p = 0.9
	#}
	s = m * (1 - p0) * p /(1-p);
	rt=c(p,s,p0);

	for(i in 1:100){
		op = p;
		os = s;
		p0=p^s;
		print(p0)
		p=m/(m2-m^2*(1-p0));
		#if (p<0.1){
		#	p = 0.1
		#}
		#if (p>=0.9){
		#	p = 0.9
		#}
		s=m * (1 - p0) * p / (1-p);
		#rt=rbind(rt,c(p,s,p0));
		rt = c(p,s,p0)
		if(abs(op-p)<0.00001 & abs(os-s)<0.00001) break;
	}
	print('change best_p0: ')
	print(p0)
	return(rt);
}

### get 0-adjusted p-value
get_pval = function(N, l, sig_0_size, sig_0_prob, num_0){
	if (N != 0){
		pval_new = pnbinom(N-1, sig_0_size, sig_0_prob, lower.tail=FALSE) / pnbinom(0, sig_0_size, sig_0_prob, lower.tail=FALSE) * (l-num_0)/l
	} else {
		pval_new = 1.0
	}
	#pval_new[pval_new>1] = 1
	return(pval_new)
}


### read data
#signal_track_file = 'MK_imm_ad.atac.meanrc.S3norm.sort.txt'
#input_track_file = 'MK_imm_ad.atac.meanrc.S3norm.sort.bg.bedgraph'
#cbg_track_file = 'atac_commonpk.txt.cbg.sort.txt'
sig = read.table(signal_track_file, header = F)
bed = sig[,1:3]
sig = sig[,4]
input = read.table(input_track_file, header = F)[,4]
cbg = scan(cbg_track_file)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
### get cbg signal track
sigcbg = sig[cbg==0]

### get sig bg regions no bgs
thresh = 0
obs_0_num = sum(sigcbg==thresh)
sigcbg_mean = mean(sigcbg)
sigcbg_moment2 = mean(sigcbg^2)
sigcbg_var = var(sigcbg)

### get zero-inflated prob & size
sigcbg_probT_sizeT = get_true_NB_prob_size(sigcbg)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sigcbg_var/sigcbg_mean, digits=3)) ))
print(sigcbg_mean)
print(sigcbg_var)
print(length(sigcbg))

### get negative binomial parameters from signal track bg regions
sigcbg_prob = sigcbg_probT_sizeT[1]
#sigcbg_prob = sigcbg_mean / sigcbg_var
if (sigcbg_prob<0.1){
	sigcbg_prob = 0.1
}

if (sigcbg_prob>=0.9){
	sigcbg_prob = 0.9
}

p0 = sigcbg_probT_sizeT[3]
#sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_mean_sig2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)
#sigcbg_size = sigcbg_mean * sigcbg_prob / (1-sigcbg_prob)
sigcbg_size = sigcbg_probT_sizeT[2]

### get input bg regions
input_mean = mean(input)
input_var = var(input)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(input_var/input_mean, digits=3)) ))
print(sigcbg_prob)
print(sigcbg_size)
print(length(input))

print(head(input))
print(summary(input))
print(input_mean)
print(input_var)

print('check NB distribution: ')
print(sigcbg_size)
print(sigcbg_prob)
print(p0)
nb_v_Test = rnbinom(1e+4, sigcbg_size, sigcbg_prob)
print(mean(nb_v_Test))
print(var(nb_v_Test))


bin_num = length(sig)[1]

### get negative binomial p-value 1st round
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) get_pval(x[1], bin_num, sigcbg_size * (x[2]+1)/(input_mean+1), sigcbg_prob, obs_0_num) )

### remove extrame p-value
nb_pval[nb_pval>1] = 1
nb_pval[nb_pval<=1e-324] = 1e-324
nb_pval_fdr = p.adjust(nb_pval, 'fdr')
nb_pval_fdr[nb_pval_fdr<=1e-324] = 1e-324
### get -log10(p-value)
neglog10_nb_pval = -log10(nb_pval)
neglog10_nb_pval_fdr = -log10(nb_pval_fdr)
print(summary(neglog10_nb_pval))
neglog10_nb_pval[neglog10_nb_pval>324] = 324
neglog10_nb_pval_fdr[neglog10_nb_pval_fdr>324] = 324
### get info
info_matrix = cbind(sigcbg_mean, sigcbg_var, sigcbg_size, sigcbg_prob)

### write output
write.table(cbind(bed, neglog10_nb_pval), paste(output_name, '.NBP.bedgraph', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
write.table(cbind(bed, neglog10_nb_pval_fdr), paste(output_name, '.NBQ.bedgraph', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
write.table(info_matrix, paste(output_name, '.mvsp.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')




