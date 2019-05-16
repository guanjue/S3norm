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

get_true_NB_prob_size_pre = function(mean_non0, mean_x2_non0){
	###### identify p0
	best_p0 = 0
	best_prob_dif = 1
	k=0
	for (i in seq(0,0.99,0.005)){
		k = k+1
		ProbT = mean_non0 / (mean_x2_non0 - mean_non0^2 * (1-i))
		if (ProbT<0.01){
			ProbT = 0.01
		}

		if (ProbT>=0.99){
			ProbT = 0.99
		}
		SizeT = mean_non0 * (1-i) * ProbT / (1-ProbT)

		nb_v_T = rnbinom(1e+4, SizeT, ProbT)
		p0_new = sum(nb_v_T==0) / length(nb_v_T)
		p0_dif = abs(i-p0_new)
		if ((k%%100)==0){
			print(paste('iteration:', toString(k)))
		}
		if (abs(i-p0_new) < best_prob_dif){
			print(paste('iteration:', toString(k)))
			print('change best_p0')
			best_prob_dif = abs(i-p0_new)
			best_p0 = i
		}
	}
	print('estimated p0: ')
	print(best_p0)
	ProbT = mean_non0 / (mean_x2_non0 - mean_non0^2 * (1-best_p0))
	SizeT = mean_non0 * (1-best_p0) * ProbT / (1-ProbT)
	return(c(ProbT, SizeT, best_p0))
}


get_true_NB_prob_size = function(x){
	m=mean(x[x>0]);
	m2=mean(x[x>0]^2);
	p0 = length(which(x==0)) / length(x);
	p = m/(m2-m^2 * (1-p0));
	if (p<0.01){
		p = 0.01
	}
	if (p>=0.99){
		p = 0.99
	}
	s = m * (1 - p0) * p /(1-p);
	rt=c(p,s,p0);

	for(i in 1:100){
		op = p;
		os = s;
		p0=p^s;
		print(p0)
		p=m/(m2-m^2*(1-p0));
		if (p<0.01){
			p = 0.01
		}
		if (p>=0.99){
			p = 0.99
		}
		s=m * (1 - p0) * p / (1-p);
		#rt=rbind(rt,c(p,s,p0));
		rt = c(p,s,p0)
		if(abs(op-p)<0.00001 & abs(os-s)<0.00001) break;
	}
	print('change best_p0: ')
	print(p0)
	return(rt);
}


get_pval_pre = function(N, l, sig_0_size, sig_0_prob, num_0){
	#3 get the probability of region that 0 < sig <= N
	pval_non0_N = pnbinom(N, sig_0_size, sig_0_prob, lower.tail = TRUE) - pnbinom(0, sig_0_size, sig_0_prob, lower.tail = TRUE)

	#4 get the probability of region that N < sig
	pval_N_ = pnbinom(N, sig_0_size, sig_0_prob, lower.tail = FALSE)

	###5 calculate the theoritical number of bins that 0 < sig <= N
	num_non0_N = l * pval_non0_N

	#6 calculate the theoritical number of bins that N < sig
	num_N_ = l * pval_N_

	#8 the new p-value for the bin
	pval_new = num_N_ / (num_0 + num_non0_N + num_N_)

	return(pval_new)
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
### make sure the min(positive number is 1)
#sig_0 = sig_0 / min(sig_0[sig_0>0])
#sig_0 = sig_0[sig_0>thresh]
obs_0_num = sum(sig_0==thresh)
sig_0_mean = mean(sig_0)
sig_0_moment2 = mean(sig_0^2)
sig_0_var = var(sig_0)


sig_0_probT_sizeT = get_true_NB_prob_size(sig_0)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_0_var/sig_0_mean, digits=3)) ))
print(sig_0_mean)
print(sig_0_var)
print(length(sig_0))

### get negative binomial parameters from signal track bg regions
sig_0_prob = sig_0_probT_sizeT[1]
#sig_0_prob = sig_0_mean / sig_0_var
if (sig_0_prob<0.1){
	sig_0_prob = 0.1
}

if (sig_0_prob>=0.9){
	sig_0_prob = 0.9
}

p0 = sig_0_probT_sizeT[3]
#sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_mean_sig2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)
#sig_0_size = sig_0_mean * sig_0_prob / (1-sig_0_prob)
sig_0_size = sig_0_probT_sizeT[2]

### get input bg regions
input_0 = input[,1]
#input_bg_non0 = input_bg[input_bg>thesh]
input_0_mean = mean(input_0)
input_0_var = var(input_0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(input_0_var/input_0_mean, digits=3)) ))
print(sig_0_prob)
print(sig_0_size)
print(length(input_0))

print(head(input_0))
print(summary(input_0))
print(input_0_mean)
print(input_0_var)

print('check NB distribution: ')
print(sig_0_size)
print(sig_0_prob)
print(p0)
nb_v_Test = rnbinom(1e+4, sig_0_size, sig_0_prob)
print(mean(nb_v_Test))
print(var(nb_v_Test))


bin_num = dim(sig)[1]

### get negative binomial p-value 1st round
nb_pval = apply(sig, MARGIN=1, function(x) get_pval(x[1], bin_num, sig_0_size, sig_0_prob, obs_0_num) )

### get -log10(p-value)
nb_pval[nb_pval<=1e-324] = 1e-324
nb_pval[nb_pval>1] = 1

### get -log10(p-value)
print('get -log10(p-value)')
print(min(nb_pval[nb_pval!=0]))
print(length(nb_pval))
print(length(nb_pval[nb_pval<=1e-324]))

### remove extrame p-value
nb_pval_min = min(nb_pval[nb_pval!=0])
print(summary(nb_pval))


############### 2nd round
thresh = 0

### get sig bg regions
sig_bg = sig[nb_pval>=0.001,]
print('sum(nb_pval>=0.001): ')
print(sum(nb_pval>=0.001))
sig_bg_non0 = sig_bg[sig_bg>thresh]
sig_bg_mean = mean(sig_bg_non0)
sig_bg_moment2 = mean(sig_bg_non0^2)
sig_bg_var = var(sig_bg_non0)

print('observed p0: ')
print(sum(sig_bg>thresh) / length(sig_bg)[1])

probT_sizeT = get_true_NB_prob_size(sig_bg)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_bg_var/sig_bg_mean, digits=3)) ))
print(sig_bg_mean)
print(sig_bg_var)
print(length(sig_bg_non0))

### get negative binomial parameters from signal track bg regions non0 regions
sig_bg_prob = probT_sizeT[1]
if (sig_bg_prob<0.01){
	sig_bg_prob = 0.01
}

if (sig_bg_prob>=0.99){
	sig_bg_prob = 0.99
}

### get p0 & size
p0 = probT_sizeT[3]
#sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_moment2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)
sig_bg_size = probT_sizeT[2]

mean_vec[1] = sig_bg_mean
var_vec[1] = sig_bg_var
size_vec[1] = sig_bg_size
prob_vec[1] = sig_bg_prob

### get negative binomial p-value

print('check NB distribution (BG): ')
print(sig_bg_size)
print(sig_bg_prob)
print(p0)
nb_v_Test = rnbinom(1e+4, sig_bg_size, sig_bg_prob)
print(mean(nb_v_Test))
print(var(nb_v_Test))


sig_input = cbind(sig, input)
#nb_pval = apply(sig_input, MARGIN=1, function(x) pnbinom(x[1], sig_bg_size * (x[2]+1)/(input_0_mean+1), sig_bg_prob, lower.tail=FALSE) )
nb_pval = apply(sig_input, MARGIN=1, function(x) get_pval(x[1], bin_num, sig_bg_size * (x[2]+1)/(input_0_mean+1), sig_bg_prob, obs_0_num) )

### get -log10(p-value)
nb_pval[nb_pval<=1e-324] = 1e-324
nb_pval[nb_pval>1] = 1
neglog10_nb_pval = -log10(nb_pval)

### write output
write.table(neglog10_nb_pval, paste(output_name, '.nbp_2r_bgadj.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

info_matrix = cbind(mean_vec, var_vec, size_vec, prob_vec)
write.table(info_matrix, paste(output_name, '.mvsp.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

