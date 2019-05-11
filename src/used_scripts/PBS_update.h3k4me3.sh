#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
#PBS -t 1-19

module load gcc/5.3.1 
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/scratch/vision/all_final_data/rc_norm
### get common pk and bg
#time Rscript get_common_pk.R rc_list.h3k4me3.txt h3k4me3_commonpk.txt 0.1
#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed h3k4me3_commonpk.txt.cpk.txt h3k4me3_commonpk.txt.cbg.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5}' | sort -k1,1 -k2,2n | cut -f4,5 > h3k4me3_commonpk.txt.cpkcbg.sort.txt
#cut -f1 h3k4me3_commonpk.txt.cpkcbg.sort.txt > h3k4me3_commonpk.txt.cpk.sort.txt
#cut -f2 h3k4me3_commonpk.txt.cpkcbg.sort.txt > h3k4me3_commonpk.txt.cbg.sort.txt

#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed h3k4me3_commonpk.txt.cpk.txt | awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1,$2,$3}' > h3k4me3_commonpk.txt.cpk.txt.bed
#time bedtools merge -i h3k4me3_commonpk.txt.cpk.txt.bed > h3k4me3_commonpk.txt.cpk.txt.merge.bed
#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed h3k4me3_commonpk.txt.cbg.txt | awk -F '\t' -v OFS='\t' '{if ($4==0) print $1,$2,$3}' > h3k4me3_commonpk.txt.cbg.txt.bed
#time bedtools merge -i h3k4me3_commonpk.txt.cbg.txt.bed > h3k4me3_commonpk.txt.cbg.txt.merge.bed

#while read LINE
#do

echo ${PBS_ARRAYID}
LINE=$(head -${PBS_ARRAYID} rc_list.h3k4me3.txt | tail -1)

sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
### get replicate mean
#paste ~/group/projects/vision/raw_5end/$sig1 ~/group/projects/vision/raw_5end/$sig2 | awk -F '\t' -v OFS='\t' '{print $1/2+$2/2}' > $ct'.'$mk'.meanrc.txt'
### get raw
cp $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.raw.txt'
### get TSnorm
time Rscript TSnorm_rc.R 'ERY_fl.h3k4me3.meanrc.txt' $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.TSnorm.txt'
### get MAnorm
time Rscript MAnorm_rc.R 'ERY_fl.h3k4me3.meanrc.txt' $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.MAnorm.txt' h3k4me3_commonpk.txt.cpk.txt
### get QTnorm
time Rscript QTnorm_rc.R 'ERY_fl.h3k4me3.meanrc.txt' $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.QTnorm.txt'
### get POIS
#time Rscript POIS_rc.R 'ERY_fl.h3k4me3.meanrc.txt' $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.POIS.txt'
### get Z
#time Rscript Z_rc.R 'ERY_fl.h3k4me3.meanrc.txt' $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.Z.txt'
### get S3norm
time python ~/group/software/S3norm/src/s3norm_cpkbg.py -r 'ERY_fl.h3k4me3.meanrc.txt' -t $ct'.'$mk'.meanrc.txt' -m 1 -i 2 -f 0.1 -n 10000 -l 10000 -a 1000 -b 0 -s /storage/home/gzx103/group/software/S3norm/src/ -p z -c F -k h3k4me3_commonpk.txt.cpk.txt -g h3k4me3_commonpk.txt.cbg.txt
cp $ct'.'$mk'.s3norm.txt' $ct'.'$mk'.meanrc.S3norm.txt'
mv $ct'.'$mk'.s3norm.txt' $ct'.'$mk'.s3norm_rc.txt'
mv $ct'.'$mk'.s3norm.scatterplot.png' $ct'.'$mk'.s3norm_rc.scatterplot.png'
mv $ct'.'$mk'.scatterplot.png' $ct'.'$mk'.rc.scatterplot.png'
mv $ct'.'$mk'.info.txt' $ct'.'$mk'.info.rc.txt'
### get NBP
# get nbp
#time Rscript ~/group/software/S3norm/src/negative_binomial_p_2r_bgadj_bayes.R $sig1 /storage/home/g/gzx103/group/projects/vision/raw_5end/ 200_noblack.5bins_bgsig_mean.round.txt /storage/home/gzx103/group/projects/vision/merged_input/ $sig1
#time Rscript ~/group/software/S3norm/src/negative_binomial_p_2r_bgadj_bayes.R $sig2 /storage/home/g/gzx103/group/projects/vision/raw_5end/ 200_noblack.5bins_bgsig_mean.round.txt /storage/home/gzx103/group/projects/vision/merged_input/ $sig2
# Fisher's method merge p
cm=$ct'.'$mk
#time Rscript ~/group/software/S3norm/src/fisher_pval.R $cm '.nbp_2r_bgadj.txt' /storage/home/gzx103/scratch/vision/all_final_data/rc_norm/ 100
### get S3norm NBP
time python ~/group/software/S3norm/src/s3norm.py -r 'ERY_fl.h3k4me3.fisher_p.txt' -t $ct'.'$mk'.fisher_p.txt' -m 1 -i 2 -f 0.1 -n 10000 -l 10000 -a 100 -b 0 -s /storage/home/gzx103/group/software/S3norm/src/ -p p -c F

### get TSnorm & raw bw
time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $ct'.'$mk'.meanrc.txt' $ct'.'$mk'.meanrc.TSnorm.txt' $ct'.'$mk'.meanrc.MAnorm.txt' $ct'.'$mk'.meanrc.QTnorm.txt' $ct'.'$mk'.meanrc.POIS.txt' $ct'.'$mk'.meanrc.Z.txt' $ct'.'$mk'.s3norm_rc.txt' $ct'.'$mk'.fisher_p.txt' $ct'.'$mk'.s3norm.txt' | sort -k1,1 -k2,2n > $ct'.'$mk'.meanrc.raw.Xnorm.txt'
time cut -f1,2,3,4 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.raw.sort.txt'
time cut -f1,2,3,5 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.TSnorm.sort.txt'
time cut -f1,2,3,6 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.MAnorm.sort.txt'
time cut -f1,2,3,7 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.QTnorm.sort.txt'
time cut -f1,2,3,8 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.POIS.sort.txt'
time cut -f1,2,3,9 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.Z.sort.txt'
time cut -f1,2,3,10 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.S3norm.sort.txt'
time cut -f1,2,3,11 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.NBP.sort.txt'
time cut -f1,2,3,12 $ct'.'$mk'.meanrc.raw.Xnorm.txt' > $ct'.'$mk'.meanrc.S3norm_NBP.sort.txt'
### get bw
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.raw.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.raw.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.TSnorm.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.TSnorm.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.MAnorm.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.MAnorm.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.QTnorm.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.QTnorm.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.POIS.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.POIS.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.Z.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.Z.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.S3norm.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.S3norm.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.NBP.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.NBP.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.S3norm_NBP.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.S3norm_NBP.bw'

#done < rc_list.h3k4me3.txt


