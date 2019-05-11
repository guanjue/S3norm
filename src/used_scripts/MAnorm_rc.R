#library(robustreg)
library(MASS)
library(affy)
#library(R.basic) # source("http://www.braju.com/R/hbLite.R") 
		 # hbLite("R.basic")
cat('test\n')
#common_peak_count_read1=read.table(input_sig1,header=FALSE)
#common_peak_count_read2=read.table(input_sig2,header=FALSE)
### get parameters
args = commandArgs(trailingOnly=TRUE)
input_sig1 = args[1]
input_sig2 = args[2]
output = args[3]
cpk_file = args[4]

small_num = 1
random_sample_num = 1000000
upperlim = 2000
lowerlim = 0

sig1 = scan(input_sig1)
sig2 = scan(input_sig2)
totalmean_sf = sum(sig1) / sum(sig2)
sig4 = (sig2+small_num) * totalmean_sf #- small_num
sig4[sig4>upperlim] = upperlim
sig4[sig4<lowerlim] = lowerlim

#sig1_notop = sig1[sig1>-1]
#sig2_notop = sig2[sig2>-1]

#sig1_z = ( sig1 - mean(sig1_notop) ) / sd(sig1_notop)
#sig2_z = ( sig1 - mean(sig2_notop) ) / sd(sig2_notop)

#sig1_z_p = pnorm(sig1_z, lower.tail = F)
#sig2_z_p = pnorm(sig2_z, lower.tail = F)

#sig1_z_p_fdr = p.adjust(sig1_z_p, 'fdr')
#sig2_z_p_fdr = p.adjust(sig2_z_p, 'fdr')

#sig1_binary = sig1_z_p_fdr < 0.05
#sig2_binary = sig2_z_p_fdr < 0.05

#sig1_binary = sig1 >= quantile(sig1, 0.99)
#sig2_binary = sig2 >= quantile(sig2, 0.99)

#peak_binary = as.logical(sig1_binary * sig2_binary) 
peak_binary = scan(cpk_file) == 1
print(sum(peak_binary))
#print(sum(sig1_binary))
#print(sum(sig2_binary))
#bg_binary_bg = as.logical((sig1_binary + sig2_binary)==0)


common_peak_count_read1 = sig1[peak_binary]+small_num
common_peak_count_read2 = sig2[peak_binary]+small_num


M=log2((common_peak_count_read2+small_num)/(common_peak_count_read1+small_num))
A=0.5*log2((common_peak_count_read2+small_num)*(common_peak_count_read1+small_num))
M = as.matrix(M)
A = as.matrix(A)

linear=lm(M~A)$coefficients
b=rlm(M~A)$coefficients


cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 = log2(common_peak_count_read1 + small_num)
log2_peak_count_read2 = log2(common_peak_count_read2 + small_num)
log2_peak_count_read2_rescaled = (2-b[2])*log2_peak_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
M_rescaled = (log2_peak_count_read2_rescaled - log2_peak_count_read1);
A_rescaled = (log2_peak_count_read2_rescaled + log2_peak_count_read1)/2;

ylim = max(c(abs(min(M)), abs(max(M)), abs(min(M_rescaled)), abs(max(M_rescaled))))

png(paste(output,".MAplot_before_rescaling.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
plot(A,M,main="MA plot before rescaling (common peaks)", ylim=c(-ylim, ylim), pch=16, cex=1)
abline(h=0,col="red",lwd=3)
abline(b,col="dodgerblue",lwd=3, lty=2)
dev.off()


png(paste(output,".MAplot_after_rescaling.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
plot(as.matrix(A_rescaled),as.matrix(M_rescaled),main=" MA plot after rescaling (all peaks)", ylim=c(-ylim, ylim), pch=16, cex=1)
abline(h=0,col="red",lwd=3)
abline(h=0,col="dodgerblue",lwd=3, lty=2)
dev.off()


#log2_allregion_count_read1 = log2(sig1 + small_num)
#log2_allregion_count_read2 = log2(sig2 + small_num)
#log2_allregion_count_read2_rescaled = (2-b[2])*log2_allregion_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
#sig2_rescaled = 2^log2_allregion_count_read2_rescaled - small_num
print(summary(sig2))
print(b)
print(((2-b[2])/(2+b[2])))
print((2*b[1]/(2+b[2])))

log2_all_count_read2_rescaled = (2-b[2])*log2(sig2+small_num)/(2+b[2]) - 2*b[1]/(2+b[2])
sig2_rescaled = 2^log2_all_count_read2_rescaled


print(summary(sig2_rescaled))


sig2_rescaled[sig2_rescaled > upperlim] = upperlim
sig2_rescaled[sig2_rescaled < lowerlim] = lowerlim


write.table(sig2_rescaled, output,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)




