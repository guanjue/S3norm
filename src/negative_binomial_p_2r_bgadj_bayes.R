library(LSD)
### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
signal_folder = args[2]
input_track_file = args[3]
input_folder = args[4]
output_name = args[5]

mean_vec = c()
var_vec = c()
size_vec = c()
prob_vec = c()


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
sig = read.table(paste(signal_folder, signal_track_file, sep=''), header = F)
input = read.table(paste(input_folder, input_track_file, sep=''), header = F)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
### get sig bg regions no bgs
thresh = 0

sig_0 = sig[,1]
sig_0_notop = sig_0[sig_0<quantile(sig_0, 0.95)]
### make sure the min(positive number is 1)
#sig_0 = sig_0 / min(sig_0[sig_0>0])
#sig_0 = sig_0[sig_0>thresh]
obs_0_num = sum(sig_0_notop==thresh)
sig_0_notop_mean = mean(sig_0_notop)
sig_0_notop_moment2 = mean(sig_0_notop^2)
sig_0_notop_var = var(sig_0_notop)


sig_0_notop_probT_sizeT = get_true_NB_prob_size(sig_0_notop)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_0_notop_var/sig_0_notop_mean, digits=3)) ))
print(sig_0_notop_mean)
print(sig_0_notop_var)
print(length(sig_0))

### get negative binomial parameters from signal track bg regions
sig_0_notop_prob = sig_0_notop_probT_sizeT[1]
#sig_0_prob = sig_0_mean / sig_0_var
if (sig_0_notop_prob<0.01){
	sig_0_notop_prob = 0.01
}

if (sig_0_notop_prob>=0.99){
	sig_0_notop_prob = 0.99
}

p0_notop = sig_0_notop_probT_sizeT[3]
#sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_mean_sig2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)
#sig_0_size = sig_0_mean * sig_0_prob / (1-sig_0_prob)
sig_0_notop_size = sig_0_notop_probT_sizeT[2]

### get input bg regions
input_0 = input[,1]
#input_bg_non0 = input_bg[input_bg>thesh]
input_0_mean = mean(input_0)
input_0_var = var(input_0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(input_0_var/input_0_mean, digits=3)) ))

print(head(input_0))
print(summary(input_0))
print(input_0_mean)
print(input_0_var)

print('check NB distribution: ')
print(sig_0_notop_size)
print(sig_0_notop_prob)
print(p0_notop)
nb_v_Test = rnbinom(1e+4, sig_0_notop_size, sig_0_notop_prob)
print(mean(nb_v_Test))
print(var(nb_v_Test))


bin_num = dim(sig)[1]

### get negative binomial p-value 1st round
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) get_pval(x[1], bin_num, sig_0_notop_size * (x[2]+1)/(input_0_mean+1), sig_0_notop_prob, obs_0_num) )

### get -log10(p-value)
nb_pval[nb_pval<=1e-323] = 1e-323
nb_pval[nb_pval>1] = 1
neglog10_nb_pval = -log10(nb_pval)

### get -log10(p-value)
print('get -log10(p-value)')
print(min(nb_pval[nb_pval!=0]))
print(length(nb_pval))
print(length(nb_pval[nb_pval<=1e-323]))

### remove extrame p-value
nb_pval_min = min(nb_pval[nb_pval!=0])
print(summary(nb_pval))


### write output
write.table(neglog10_nb_pval, paste(output_name, '.nbp_2r_bgadj.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
mean_vec[1] = sig_0_notop_mean
var_vec[1] = sig_0_notop_var
size_vec[1] = sig_0_notop_size
prob_vec[1] = sig_0_notop_prob

info_matrix = cbind(mean_vec, var_vec, size_vec, prob_vec)
write.table(info_matrix, paste(output_name, '.mvsp.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


