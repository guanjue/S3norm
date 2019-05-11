#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
#PBS -t 1-18

module load gcc/5.3.1 
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/scratch/vision/all_final_data/rc_norm
#while read LINE
#do

echo ${PBS_ARRAYID}
LINE=$(head -${PBS_ARRAYID} rc_list.ctcf.txt | tail -1)

sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
cm=$ct'.'$mk

### get rc_norm matrix
time cut -f4,5,6,7,10 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.raw.Xnorm_rc.txt'
time Rscript ~/scratch/vision/all_final_data/rc_norm_rep/get_frip.R $ct'.'$mk'.meanrc.raw.Xnorm_rc.txt' $ct'.'$mk'.frip.txt' 
time Rscript ~/scratch/vision/all_final_data/rc_norm_rep/get_mean_sig_bgpk.R $ct'.'$mk'.meanrc.raw.Xnorm_rc.txt' $ct'.'$mk'.pkbg.txt'


#done < rc_list.ctcf.txt


