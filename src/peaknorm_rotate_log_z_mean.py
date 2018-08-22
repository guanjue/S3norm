#module load python/2.7
import os
import numpy as np
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

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
		print "Method is not implemented"
		p0 = None
	return p0

################################################################################################
### NewtonRaphsonMethod
def NewtonRaphsonMethod(sig1_pk,sig1_bg, sig2_pk,sig2_bg, A,B, moment, converge_thresh, numIterations):
	sig1_pk_mean = np.mean(sig1_pk**moment)
	sig1_bg_mean = np.mean(sig1_bg**moment)

	for i in range(0, numIterations):
		fb = sig1_bg_mean * np.mean(sig2_pk**(moment*B)) - sig1_pk_mean * np.mean(sig2_bg**(moment*B))
		dfb = moment * sig1_bg_mean * np.mean(np.log(sig2_pk) * sig2_pk**(moment*B)) - moment * sig1_pk_mean * np.mean(np.log(sig2_bg) * sig2_bg**(moment*B))

		### next step
		B = B - fb / dfb	
		A = sig1_bg_mean / np.mean(sig2_bg**(moment*B))

		print("Iteration %d | dFB: %f" % (i, dfb))
		print([A,B])

		last_AB = [A, B]

		if abs(fb / dfb) < converge_thresh:
			print('converged!')
			used_AB = [A, B]
			print('pk')
			print(sig1_pk_mean)
			print(np.mean(sig2_bg))
			print('bg')
			print(sig1_bg_mean)
			print(np.mean(sig2_bg))
			break

	if abs(fb / dfb) >= converge_thresh:
		print('NOT converged...')
		used_AB = last_AB

	print('used')
	print(used_AB)
	return np.array(used_AB)

################################################################################################
###
def pknorm(sample_num, sig1_wg_raw, sig2_wg_raw, upperlim, lowerlim):
	sig1_output_name = sig1_wg_raw.split('.')[0]+'.'+sig1_wg_raw.split('.')[1]
	sig2_output_name = sig2_wg_raw.split('.')[0]+'.'+sig2_wg_raw.split('.')[1]

	### add small_number
	small_num = 1e-1

	### read whole genome signals
	sig1 = read2d_array(sig1_wg_raw, float)
	sig2 = read2d_array(sig2_wg_raw, float)

	### total reads norm
	if sig1_output_name == sig2_output_name:
		sig2 = sig1
	
	### read whole genome binary label
	#sig1_z_p_fdr = p_adjust(1 - norm.cdf((sig1 - np.mean(sig1))/ np.std(sig1)), 'fdr')
	#sig1_binary = sig1_z_p_fdr < 0.05
	#sig1_pk_num = np.sum(sig1_binary)
	#if sig1_pk_num <= 1e4:
	#	sig1_thresh = np.sort(sig1, axis=None)[-10000]
	#	print('rank sig1')
	#	sig1_binary = sig1 > sig1_thresh
	sig1_binary = 10**(-sig1) <= 0.001
	#print(sig1_pk_num)

	#sig2_z_p_fdr = p_adjust(1 - norm.cdf((sig2 - np.mean(sig2))/ np.std(sig2)), 'fdr')
	#sig2_binary = sig2_z_p_fdr < 0.05
	#sig2_pk_num = np.sum(sig2_binary)
	#if sig2_pk_num <= 1e4:
	#	sig2_thresh = np.sort(sig2, axis=None)[-10000]
	#	print('rank sig2')
	#	sig2_binary = sig2 > sig2_thresh
	sig2_binary = 10**(-sig2) <= 0.001
	#print(sig2_pk_num)

	### peak region (both != 0 in sig1 & sig2)
	peak_binary_pk = (sig1_binary[:,0] & sig2_binary[:,0])
	print(np.sum(peak_binary_pk))
	peak_binary = peak_binary_pk & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) 
	print(np.max(sig1[:,0]))
	print(np.sum(peak_binary))

	### background region (both == 0 in sig1 & sig2)
	bg_binary_bg = ~(sig1_binary[:,0] | sig2_binary[:,0])
	print(np.sum(bg_binary_bg))
	bg_binary = bg_binary_bg & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) 
	print(np.sum(bg_binary))


	### get transformation factor
	AB = NewtonRaphsonMethod(sig1[peak_binary,0]+small_num,sig1[bg_binary,0]+small_num, sig2[peak_binary,0]+small_num,sig2[bg_binary,0]+small_num, 1.0, 2.0, 1, 0.0001, 200)
	A=AB[0]
	B=AB[1]
	print('transformation: '+'B: '+str(B)+'; A: '+str(A))
	### transformation
	sig2_norm = []
	for s in sig2[:,0]:
		s = s
		s_norm = (A* (s+small_num)**B) - small_num
		if s_norm >= upperlim:
			s_norm = upperlim
		elif s_norm <= lowerlim:
			s_norm = lowerlim
		sig2_norm.append(s_norm)

	### total reads sf (for compare)
	sig1_totalmean = np.mean(sig1[(sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]),0])
	sig2_totalmean = np.mean(sig2[(sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]),0])
	total_mean_sf = sig1_totalmean / sig2_totalmean

	### convert to float np.array
	sig2_norm = np.array(sig2_norm, float)
	sig2_norm_totalmean = np.mean(sig2_norm[(sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0])])
	print('total means: ')
	print(sig1_totalmean)
	print(sig2_totalmean)
	print(sig2_norm_totalmean)
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

	### write output: normalized signal
	write2d_array(sig2_norm, sig2_output_name + '.pknorm.txt')

	### write output: sf & FRiP
	info = np.array([[total_mean_sf, B, A], [sig1_FRiP, sig2_norm_FRiP, sig2_FRiP]])
	write2d_array(info, sig2_output_name + '.info.txt')


	### plot scatter plot
	np.random.seed(2018)
	idx = np.random.randint(sig2_norm.shape[0], size=sample_num)
	peak_binary_sample = peak_binary_pk[idx]
	bg_binary_sample = bg_binary_bg[idx]
	plot_x = np.log2(sig2_norm[idx,0]+small_num)
	plot_y = np.log2(sig1[idx,0]+small_num)
	plot_xn = np.log2(sig2[idx,0]+small_num)
	plot_yn = np.log2(sig1[idx,0]+small_num)
	plot_totalmean_xn = np.log2(sig2_norm_totalmean+small_num)
	plot_totalmean_yn = np.log2(sig1_totalmean+small_num)

	lims_max = np.max(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))
	lims_min = np.min(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))

	plt.figure()
	plt.scatter(plot_x, plot_y, marker='.', color='dodgerblue')
	plt.scatter(plot_x[bg_binary_sample], plot_y[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.scatter(sig2_1log_pk_m_pkn, sig1_1log_pk_m_od, marker='.', color='k')
	plt.scatter(sig2_1log_bg_m_pkn, sig1_1log_bg_m_od, marker='.', color='k')
	plt.scatter(plot_totalmean_xn, plot_totalmean_yn, marker='.', color='red')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_pkn, sig2_1log_pk_m_pkn], [sig1_1log_bg_m_od, sig1_1log_pk_m_od])
	#plt.scatter(np.mean(plot_x[peak_binary_sample]), np.mean(plot_y[peak_binary_sample]), marker='.', color='k')
	#plt.scatter(np.mean(plot_x[bg_binary_sample]), np.mean(plot_y[bg_binary_sample]), marker='.', color='k')
	#plt.plot([np.mean(plot_x[bg_binary_sample]), np.mean(plot_x[peak_binary_sample])], [np.mean(plot_y[bg_binary_sample]), np.mean(plot_y[peak_binary_sample])])
	plt.xlabel(sig2_output_name + '.pknorm')
	plt.ylabel(sig1_output_name + '.pknorm')
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.axis('scaled')
	plt.savefig(sig2_output_name + '.pknorm.scatterplot.png')


	plot_xn = np.log2(sig2[idx,0]+small_num)
	plot_yn = np.log2(sig1[idx,0]+small_num)
	plot_totalmean_xn = np.log2(sig2_totalmean+small_num)
	plot_totalmean_yn = np.log2(sig1_totalmean+small_num)

	lims_max = np.max(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))
	lims_min = np.min(np.concatenate((plot_x, plot_y, plot_xn, plot_yn)))

	plt.figure()
	plt.scatter(plot_xn, plot_yn, marker='.', color='dodgerblue')
	plt.scatter(plot_xn[bg_binary_sample], plot_yn[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_xn[peak_binary_sample], plot_yn[peak_binary_sample], marker='.', color='coral')
	plt.scatter(sig2_1log_pk_m_od, sig1_1log_pk_m_od, marker='.', color='k')
	plt.scatter(sig2_1log_bg_m_od, sig1_1log_bg_m_od, marker='.', color='k')
	plt.scatter(plot_totalmean_xn, plot_totalmean_yn, marker='.', color='red')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_od, sig2_1log_pk_m_od], [sig1_1log_bg_m_od, sig1_1log_pk_m_od])
	#plt.scatter(np.mean(plot_xn[peak_binary_sample]), np.mean(plot_yn[peak_binary_sample]), marker='.', color='k')
	#plt.scatter(np.mean(plot_xn[bg_binary_sample]), np.mean(plot_yn[bg_binary_sample]), marker='.', color='k')
	#plt.plot([np.mean(plot_xn[bg_binary_sample]), np.mean(plot_xn[peak_binary_sample])], [np.mean(plot_yn[bg_binary_sample]), np.mean(plot_yn[peak_binary_sample])])
	plt.xlabel(sig2_output_name + '.pknorm')
	plt.ylabel(sig1_output_name + '.pknorm')
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.axis('scaled')
	plt.savefig(sig2_output_name + '.scatterplot.png')


############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hn:a:b:u:l:")
	except getopt.GetoptError:
		print 'time python peaknorm_rotate_log_z_mean.py -n sample_num -a reference_dataset -b target_dataset -u upperlim -l lowerlim'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python peaknorm_rotate_log_z_mean.py -n sample_num -a reference_dataset -b target_dataset -u upperlim -l lowerlim'
			sys.exit()
		elif opt=="-n":
			sample_num=int(arg.strip())		
		elif opt=="-a":
			sig1_wg_raw=str(arg.strip())					
		elif opt=="-b":
			sig2_wg_raw=str(arg.strip())		
		elif opt=="-u":
			upperlim=float(arg.strip())
		elif opt=="-l":
			lowerlim=float(arg.strip())

	pknorm(sample_num, sig1_wg_raw, sig2_wg_raw, upperlim, lowerlim)

if __name__=="__main__":
	main(sys.argv[1:])


