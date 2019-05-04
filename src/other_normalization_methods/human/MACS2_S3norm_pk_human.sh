#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
#PBS -t 1-20

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/scratch/vision/all_final_data/S3norm_pk_human

echo ${PBS_ARRAYID}

LINE=$(head -${PBS_ARRAYID} atac_file_list.txt | tail -1)
ct=$(echo "$LINE" | awk '{print $1}')
echo $ct

### get bed coordinates
#cat windowsNoBlackForIdeas.bed | awk -F ' ' -v OFS='\t' '{print $1, $2, $3}' > human_vision_200bp.bed
#cat windowsNoBlackForIdeas.bed | awk -F ' ' -v OFS='\t' '{if ($1!="chrEBV") print $1, $2, $3, $1"_"$2"_"$3}' > human_vision.200bp.bed
#cat windowsNoBlackForIdeas.bed | awk -F ' ' -v OFS='\t' '{if ($1!="chrEBV" && int(($2+$3)/2)-2500 > 0) print $1, int(($2+$3)/2)-2500, int(($2+$3)/2)+2500, $1"_"$2"_"$3; else if ($1!="chrEBV" && int(($2+$3)/2)-2500 <= 0) print $1, 0, int(($2+$3)/2)+2500, $1"_"$2"_"$3}' > human_vision.5kb.bed
#cat windowsNoBlackForIdeas.bed | awk -F ' ' -v OFS='\t' '{if ($1!="chrEBV" && int(($2+$3)/2)-5000 > 0) print $1, int(($2+$3)/2)-5000, int(($2+$3)/2)+5000, $1"_"$2"_"$3; else if ($1!="chrEBV" && int(($2+$3)/2)-5000 <= 0) print $1, 0, int(($2+$3)/2)+5000, $1"_"$2"_"$3}' > human_vision.10kb.bed

### get bedgraph
time paste human_vision_200bp.bed $ct'.ATAC.s3norm.0_100.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n > $ct'.ATAC.s3norm.0_100.bedgraph'

###
time ~/group/software/ucsc/bedGraphToBigWig $ct'.ATAC.s3norm.0_100.bedgraph' ~/group/genome/hg38/hg38.chrom.sizes $ct'.ATAC.s3norm.0_100.bw'

### bigWigAverageOverBed
#time ~/group/software/ucsc/bigWigAverageOverBed $ct'.ATAC.s3norm.0_100.bw' ~/group/projects/vision/merged_input/200_noblack.11_22_2017.5kb.bed $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab'


