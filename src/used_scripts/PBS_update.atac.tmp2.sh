#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
#PBS -t 1-18

module load gcc/5.3.1 
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/scratch/vision/all_final_data/rc_norm
### get common pk and bg
#time Rscript get_common_pk.R rc_list.atac.txt atac_commonpk.txt 0.1
#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed atac_commonpk.txt.cpk.txt atac_commonpk.txt.cbg.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5}' | sort -k1,1 -k2,2n | cut -f4,5 > atac_commonpk.txt.cpkcbg.sort.txt
#cut -f1 atac_commonpk.txt.cpkcbg.sort.txt > atac_commonpk.txt.cpk.sort.txt
#cut -f2 atac_commonpk.txt.cpkcbg.sort.txt > atac_commonpk.txt.cbg.sort.txt

#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed atac_commonpk.txt.cpk.txt | awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1,$2,$3}' > atac_commonpk.txt.cpk.txt.bed
#time bedtools merge -i atac_commonpk.txt.cpk.txt.bed > atac_commonpk.txt.cpk.txt.merge.bed
#time paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed atac_commonpk.txt.cbg.txt | awk -F '\t' -v OFS='\t' '{if ($4==0) print $1,$2,$3}' > atac_commonpk.txt.cbg.txt.bed
#time bedtools merge -i atac_commonpk.txt.cbg.txt.bed > atac_commonpk.txt.cbg.txt.merge.bed


#while read LINE
#do

echo ${PBS_ARRAYID}

LINE=$(head -${PBS_ARRAYID} rc_list.atac.txt | tail -1)

sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
cm=$ct'.'$mk
time cut -f1,2,3,4 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.raw.sort.10k.txt'
time cut -f1,2,3,5 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.TSnorm.sort.10k.txt'
time cut -f1,2,3,6 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.MAnorm.sort.10k.txt'
time cut -f1,2,3,7 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.QTnorm.sort.10k.txt'
#time cut -f1,2,3,8 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.POIS.sort.txt'
#time cut -f1,2,3,9 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.Z.sort.txt'
time cut -f1,2,3,10 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.S3norm.sort.10k.txt'
#time cut -f1,2,3,11 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.NBP.sort.txt'
#time cut -f1,2,3,12 $ct'.'$mk'.meanrc.raw.Xnorm.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=100000) print $1,$2,$3,$4; else print $1,$2,$3,100000}' > $ct'.'$mk'.meanrc.S3norm_NBP.sort.txt'
### get bw
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.raw.sort.10k.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.raw.10k.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.TSnorm.sort.10k.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.TSnorm.10k.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.MAnorm.sort.10k.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.MAnorm.10k.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.QTnorm.sort.10k.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.QTnorm.10k.bw'
#time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.POIS.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.POIS.bw'
#time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.Z.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.Z.bw'
time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.S3norm.sort.10k.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.S3norm.10k.bw'
#time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.NBP.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.NBP.bw'
#time ~/group/software/ucsc/bedGraphToBigWig $ct'.'$mk'.meanrc.S3norm_NBP.sort.txt' ~/group/projects/vision/input_norm/mm10.chrom.sizes $ct'.'$mk'.meanrc.S3norm_NBP.bw'

#done < rc_list.atac.txt


