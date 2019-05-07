#module load python/2.7
import os
from subprocess import call
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm, nbinom

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### p-value adjust (fdr & bonferroni)
def p_adjust(pvalue, method):
	p = pvalue
	n = len(p)
	p0 = np.copy(p, order='K')
	nna = np.isnan(p)
	p = p[~nna]
	lp = len(p)
	if method == "bonferroni":
		p0[~nna] = np.fmin(1, lp * p)
	elif method == "fdr":
		i = np.arange(lp, 0, -1)
		o = (np.argsort(p))[::-1]
		ro = np.argsort(o)
		p0[~nna] = np.fmin(1, np.minimum.accumulate((p[o]/i*lp)))[ro]
	else:
		print("Method is not implemented")
		p0 = None
	return p0

################################################################################################
### NewtonRaphsonMethod
def NewtonRaphsonMethod(sig1_pk,sig1_bg, sig2_pk,sig2_bg, upperlim, A,B, moment, converge_thresh, numIterations):
	np.random.seed(2018)
	sig1_pk_mean = np.mean(sig1_pk**moment)
	sig1_bg_mean = np.mean(sig1_bg**moment)
	converge_min = 1000.0
	for i in range(0, numIterations):
		sig2_pk_transformed = sig2_pk**(moment*B)
		sig2_bg_transformed = sig2_bg**(moment*B)
		sig2_pk_transformed[sig2_pk_transformed>upperlim] = upperlim
		sig2_bg_transformed[sig2_bg_transformed>upperlim] = upperlim
		fb = sig1_bg_mean * np.mean(sig2_pk_transformed) - sig1_pk_mean * np.mean(sig2_bg_transformed)
		dfb = moment * sig1_bg_mean * np.mean(np.log(sig2_pk) * sig2_pk_transformed) - moment * sig1_pk_mean * np.mean(np.log(sig2_bg) * sig2_bg_transformed)

		### next step
		B = B - fb / dfb
		if B <=0:
			B = 1.0 + np.random.normal(0,1,1)[0]
		sig2_bg_transformed = sig2_bg**(moment*B)
		sig2_bg_transformed[sig2_bg_transformed>upperlim] = upperlim
		A = sig1_bg_mean / np.mean(sig2_bg_transformed)

		print("Iteration %d | dFB: %f" % (i, fb/dfb))
		print([A,B])

		if abs(fb / dfb) < converge_min:
			i_min = i
			converge_min = abs(fb / dfb)
			converge_min_AB = [A, B]

		if abs(fb / dfb) < converge_thresh:
			print('converged!')
			used_AB = [A, B]
			break

	if abs(fb / dfb) >= converge_thresh:
		print('NOT converged...')
		print("Iteration %d | dFB: %f" % (i_min, converge_min))
		B = ( np.mean(np.log(sig1_pk)) - np.mean(np.log(sig1_bg)) ) / ( np.mean(np.log(sig2_pk)) - np.mean(np.log(sig2_bg)) )
		sig2_bg_transformed = sig2_bg**(moment*B)
		sig2_bg_transformed[sig2_bg_transformed>upperlim] = upperlim
		A = sig1_bg_mean / np.mean(sig2_bg_transformed)
		mean_log_AB = [A, B]
		used_AB = mean_log_AB
		print('converge_min_AB')
		print(converge_min_AB)

	print('used: ')
	print(used_AB)
	return np.array(used_AB)

################################################################################################
### Negative binomial p-value
def nb_cpf(signal_vec):
	sig_mean = np.mean(signal_vec)
	sig_var = np.var(signal_vec)
	sig_prob = sig_mean / sig_var
	if sig_prob < 0.1:
		sig_prob = 0.1
	elif sig_prob > 0.9:
		sig_prob = 0.9
	sig_size = sig_mean * sig_prob / (1-sig_prob)
	nbp = 1-nbinom.cdf(signal_vec, sig_size, sig_prob)
	return nbp

################################################################################################
### s3norm
def s3norm(sig1_wg_raw, sig2_wg_raw, moment, B_init, fdr_thresh, sample_num, rank_lim, upperlim, lowerlim, script_folder, p_method, ref, common_pk_binary, common_bg_binary):
	sig1_output_name = sig1_wg_raw.split('.')[0]+'.'+sig1_wg_raw.split('.')[1]
	sig2_output_name = sig2_wg_raw.split('.')[0]+'.'+sig2_wg_raw.split('.')[1]

	### read whole genome signals
	sig1 = read2d_array(sig1_wg_raw, float)
	sig2 = read2d_array(sig2_wg_raw, float)

	### limit signal
	sig1[sig1>upperlim] = upperlim
	sig1[sig1<lowerlim] = lowerlim
	sig2[sig2>upperlim] = upperlim
	sig2[sig2<lowerlim] = lowerlim

	### read whole genome binary label
	if p_method == 'z':
		#sig1_log2 = np.log2(sig1+0.01)
		#sig1_z_p_fdr = p_adjust(1 - norm.cdf((sig1_log2 - np.mean(sig1_log2))/ np.std(sig1_log2)), 'fdr')
		sig1a = sig1+0.01
		sig1_z_p_fdr = p_adjust(1 - norm.cdf((sig1a - np.mean(sig1a))/ np.std(sig1a)), 'fdr')
		sig1_binary = sig1_z_p_fdr < fdr_thresh
	elif p_method == 't':
		sig1_binary = sig1 >= np.quantile(sig1, fdr_thresh)
	elif p_method == 'p':
		sig1_p = 10.0**(-sig1)
		sig1_z_p_fdr = p_adjust(sig1_p, 'fdr')
		write2d_array(sig1_z_p_fdr.reshape(sig1_z_p_fdr.shape[0],1), sig1_wg_raw + '.nbp.txt')
		sig1_binary = (sig1_z_p_fdr < fdr_thresh)

	sig1_pk_num = np.sum(sig1_binary)
	print('sig1_pk_num')
	print(sig1_pk_num)

	### if pk number < 10000 then use rank for sig2
	if sig1_pk_num <= rank_lim:
		sig1_thresh = np.sort(sig1, axis=None)[-rank_lim]
		print('rank sig1')
		sig1_binary = sig1 > sig1_thresh

	print('sig1_pk_num')
	print(sum(sig1_binary))
	print(sig1_pk_num)

	if p_method == 'z':
		#sig2_log2 = np.log2(sig2+0.01)
		#sig2_z_p_fdr = p_adjust(1 - norm.cdf((sig2_log2 - np.mean(sig2_log2))/ np.std(sig2_log2)), 'fdr')
		sig2a = sig2+0.01
		sig2_z_p_fdr = p_adjust(1 - norm.cdf((sig2a - np.mean(sig2a))/ np.std(sig2a)), 'fdr')
		sig2_binary = sig2_z_p_fdr < fdr_thresh
	elif p_method == 't':
		sig2_binary = sig2 >= np.quantile(sig2, fdr_thresh)
	elif p_method == 'p':
		sig2_p = 10.0**(-sig2)
		sig2_z_p_fdr = p_adjust(sig2_p, 'fdr')
		write2d_array(sig2_z_p_fdr.reshape(sig2_z_p_fdr.shape[0],1), sig2_wg_raw + '.nbp.txt')
		sig2_binary = (sig2_z_p_fdr < fdr_thresh)

	sig2_pk_num = np.sum(sig2_binary)
	print('sig2_pk_num')
	print(sig2_pk_num)

	### if pk number < 10000 then use rank for sig2
	if sig2_pk_num <= rank_lim:
		sig2_thresh = np.sort(sig2, axis=None)[-rank_lim]
		print('rank sig2')
		sig2_binary = sig2 > sig2_thresh

	print('sig2_pk_num')
	print(sum(sig2_binary))
	print(sig2_pk_num)

	if (common_pk_binary!='0') and (common_bg_binary!='0'):
		common_pk_binary_sig = read2d_array(common_pk_binary, float)
		peak_binary = (common_pk_binary_sig == 1.0)[:,0]
		common_bg_binary_sig = read2d_array(common_bg_binary, float)
		bg_binary = (common_bg_binary_sig == 0.0)[:,0]

	### peak region (both != 0 in sig1 & sig2)
	if ref == 'F':
		peak_binary = (sig1_binary[:,0] & sig2_binary[:,0])
	else:
		peak_binary = (sig1_binary[:,0] | sig2_binary[:,0])
	print(np.sum(peak_binary))

	### background region (both == 0 in sig1 & sig2)
	bg_binary = ~(sig1_binary[:,0] | sig2_binary[:,0])
	print(np.sum(bg_binary))

	### get common bg pk
	sig1_cbg = sig1[bg_binary,0]
	sig2_cbg = sig2[bg_binary,0]
	sig1_cpk = sig1[peak_binary,0]
	sig2_cpk = sig2[peak_binary,0]

	### get data driven added small number
	sig1_cbg_mean = np.mean(sig1_cbg)
	sig2_cbg_mean = np.mean(sig2_cbg)
	sig1_cpk_mean = np.mean(sig1_cpk)
	sig2_cpk_mean = np.mean(sig2_cpk)

	if (sig1_output_name != sig2_output_name) and (sig1_cbg_mean != sig2_cbg_mean) and (sig1_cpk_mean != sig2_cpk_mean):
		small_num = (sig1_cpk_mean*sig2_cbg_mean - sig1_cbg_mean*sig2_cpk_mean) / ((sig1_cbg_mean-sig1_cpk_mean)-(sig2_cbg_mean-sig2_cpk_mean))
		#small_num = 1.0
	else:
		small_num = 0.0

	if small_num >1:
		small_num = 1.0
	elif small_num <0.1:
		small_num = 0.1
	print('added small number: '+str(small_num))

	### get transformation factor
	if sig1_output_name != sig2_output_name:
		AB = NewtonRaphsonMethod(sig1_cpk+small_num,sig1_cbg+small_num, sig2_cpk+small_num,sig2_cbg+small_num, upperlim, 0.5, 2.0, moment, 1e-5, 200)
		A=AB[0]
		B=AB[1]
	else:
		A=1.0
		B=1.0		
	print('transformation: '+'B: '+str(B)+'; A: '+str(A))
	### transformation
	sig2_norm = []
	for s in sig2[:,0]:
		s_norm = (A * (s+small_num)**B) - small_num
		if s_norm > upperlim:
			s_norm = upperlim
		elif s_norm < lowerlim:
			s_norm = lowerlim
		sig2_norm.append(s_norm)
	sig2_norm = np.array(sig2_norm)

	### get transformation factor R2
	sig2_norm[sig2_norm>upperlim] = upperlim
	sig2_norm[sig2_norm<lowerlim] = lowerlim
	sig2_cbg = sig2_norm[bg_binary]
	sig2_cpk = sig2_norm[peak_binary]
	if sig1_output_name != sig2_output_name:
		AB = NewtonRaphsonMethod(sig1_cpk+small_num,sig1_cbg+small_num, sig2_cpk+small_num,sig2_cbg+small_num, upperlim, A, B, moment, 1e-5, 200)
		A=AB[0]
		B=AB[1]
	else:
		A=1.0
		B=1.0
	print('R2 transformation: '+'B: '+str(B)+'; A: '+str(A))
	### transformation
	sig2_norm_R2 = []
	for s in sig2_norm:
		s_norm = (A * (s+small_num)**B) - small_num
		if s_norm > upperlim:
			s_norm = upperlim
		elif s_norm < lowerlim:
			s_norm = lowerlim
		sig2_norm_R2.append(s_norm)
	sig2_norm = np.array(sig2_norm_R2)


	### total reads sf (for compare)
	sig1_totalmean = np.mean(sig1)
	sig2_totalmean = np.mean(sig2)
	total_mean_sf = sig1_totalmean / sig2_totalmean

	### convert to float np.array
	sig2_norm = np.array(sig2_norm, float)
	sig2_norm[sig2[:,0]==0] = 0.0
	sig2_norm_totalmean = np.mean(sig2_norm)
	### reshape for writing oputput
	sig2_norm = np.reshape(sig2_norm, (sig2_norm.shape[0],1))

	### rotated means for sig2 for plotting
	sig1_1log_pk_m_od = np.log2(np.mean(sig1[peak_binary,0])+small_num)
	sig2_1log_pk_m_od = np.log2(np.mean(sig2[peak_binary,0])+small_num)

	sig1_1log_bg_m_od = np.log2(np.mean(sig1[bg_binary,0])+small_num)
	sig2_1log_bg_m_od = np.log2(np.mean(sig2[bg_binary,0])+small_num)

	sig2_1log_pk_m_pkn = np.log2(np.mean(sig2_norm[peak_binary,0])+small_num)
	sig2_1log_bg_m_pkn = np.log2(np.mean(sig2_norm[bg_binary,0])+small_num)

	###FRiP score
	sig2_norm_FRiP = np.sum(sig2_norm[(sig2_binary[:,0]!=0),0]) / np.sum(sig2_norm)
	sig2_FRiP = np.sum(sig2[(sig2_binary[:,0]!=0),0]) / np.sum(sig2)
	sig1_FRiP = np.sum(sig1[(sig1_binary[:,0]!=0),0]) / np.sum(sig1)

	if ref == 'F':
		### write output: normalized signal
		write2d_array(sig2_norm, sig2_output_name + '.s3norm.txt')
		### write output: sf & FRiP
		info = np.array([[total_mean_sf, B, A], [sig1_FRiP, sig2_norm_FRiP, sig2_FRiP]])
		write2d_array(info, sig2_output_name + '.info.txt')
	else:
		### write output: normalized signal
		write2d_array(sig2_norm, sig2_output_name + '.s3norm.ref.txt')
		### write output: sf & FRiP
		info = np.array([[total_mean_sf, B, A], [sig1_FRiP, sig2_norm_FRiP, sig2_FRiP]])
		write2d_array(info, sig2_output_name + '.ref.info.txt')


	### plot scatter plot
	np.random.seed(2018)
	idx = np.random.randint(sig2_norm.shape[0], size=sample_num)
	peak_binary_sample = peak_binary[idx]
	bg_binary_sample = bg_binary[idx]
	plot_x = np.log2(sig2_norm[idx,0]+small_num)
	plot_y = np.log2(sig1[idx,0]+small_num)
	plot_xn = np.log2(sig2[idx,0]+small_num)
	plot_yn = np.log2(sig1[idx,0]+small_num)
	lims_max = np.max(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))
	lims_min = np.min(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))

	plt.figure(figsize=(5, 5))
	plt.scatter(plot_x, plot_y, marker='.', color='dodgerblue')
	plt.scatter(plot_x[bg_binary_sample], plot_y[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_pkn, sig2_1log_pk_m_pkn], [sig1_1log_bg_m_od, sig1_1log_pk_m_od], linestyle='dashed', color='lime')
	plt.scatter(sig2_1log_pk_m_pkn, sig1_1log_pk_m_od, marker='.', color='red', s=60)
	plt.scatter(sig2_1log_bg_m_pkn, sig1_1log_bg_m_od, marker='.', color='black', s=60)
	plt.scatter(np.log2(sig2_norm_totalmean), np.log2(sig1_totalmean), marker='.', color='blue', s=60)
	plt.xlabel(sig2_output_name + '.s3norm', fontsize=12)
	plt.ylabel(sig1_output_name + '.s3norm', fontsize=12)
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.axis('scaled')
	plt.savefig(sig2_output_name + '.s3norm.scatterplot.png')


	plot_xn = np.log2(sig2[idx,0]+small_num)
	plot_yn = np.log2(sig1[idx,0]+small_num)
	lims_max = np.max(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))
	lims_min = np.min(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))

	plt.figure(figsize=(5, 5))
	plt.scatter(plot_xn, plot_yn, marker='.', color='dodgerblue')
	plt.scatter(plot_xn[bg_binary_sample], plot_yn[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_xn[peak_binary_sample], plot_yn[peak_binary_sample], marker='.', color='coral')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_od, sig2_1log_pk_m_od], [sig1_1log_bg_m_od, sig1_1log_pk_m_od], linestyle='dashed', color='lime')
	plt.scatter(sig2_1log_pk_m_od, sig1_1log_pk_m_od, marker='.', color='red', s=60)
	plt.scatter(sig2_1log_bg_m_od, sig1_1log_bg_m_od, marker='.', color='black', s=60)
	plt.scatter(np.log2(sig2_totalmean), np.log2(sig1_totalmean), marker='.', color='blue', s=60)
	plt.xlabel(sig2_output_name + '.s3norm', fontsize=12)
	plt.ylabel(sig1_output_name + '.s3norm', fontsize=12)
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.axis('scaled')
	plt.savefig(sig2_output_name + '.scatterplot.png')


############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hr:t:m:i:f:n:l:a:b:s:p:c:k:g:")
	except getopt.GetoptError:
		print('time python s3norm.py -r reference_signal_track.txt -t target_signal_track.txt -m moment -i initial_B -f fdrthresh -n plotpoints_num -l rank_lim -a upperlimit -b lowerlimit -s script_folder-p p-value_method -c used_ref_version -k common_pk_binary -g common_bg_binary')
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print('time python s3norm.py -r reference_signal_track.txt -t target_signal_track.txt -m moment -i initial_B -f fdrthresh -n plotpoints_num -l rank_lim -a upperlimit -b lowerlimit -s script_folder -p p-value_method -c used_ref_version -k common_pk_binary -g common_bg_binary')	
		elif opt=="-r":
			sig1_wg_raw=str(arg.strip())				
		elif opt=="-t":
			sig2_wg_raw=str(arg.strip())
		elif opt=="-m":
			moment=int(arg.strip())
		elif opt=="-i":
			B_init=float(arg.strip())
		elif opt=="-f":
			fdr_thresh=float(arg.strip())
		elif opt=="-n":
			sample_num=int(arg.strip())
		elif opt=="-l":
			rank_lim=int(arg.strip())
		elif opt=="-a":
			upperlim=float(arg.strip())
		elif opt=="-b":
			lowerlim=float(arg.strip())
		elif opt=="-s":
			script_folder=str(arg.strip())
		elif opt=="-p":
			p_method=str(arg.strip())
		elif opt=="-c":
			ref=str(arg.strip())
		elif opt=="-k":
			common_pk_binary=str(arg.strip())
		elif opt=="-g":
			common_bg_binary=str(arg.strip())

	s3norm(sig1_wg_raw, sig2_wg_raw, moment, B_init, fdr_thresh, sample_num, rank_lim, upperlim, lowerlim, script_folder, p_method, ref, common_pk_binary, common_bg_binary)

if __name__=="__main__":
	main(sys.argv[1:])


