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

echo ${PBS_ARRAYID}

#declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "POIS" "Z" "S3norm" "NBP" "S3norm_NBP")
declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

LINE=$(head -${PBS_ARRAYID} rc_list.atac.txt | tail -1)
sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
mkdir $ct'_'$mk'_macs_bg'
### get pk calling result and pk number
for type in "${type_vector[@]}"
do
	echo $type
	### get bg signals
	time ~/group/software/ucsc/bigWigAverageOverBed $ct'.'$mk'.meanrc.'$type'.bw' ~/group/projects/vision/merged_input/200_noblack.11_22_2017.5kb.bed $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab'
	time ~/group/software/ucsc/bigWigAverageOverBed $ct'.'$mk'.meanrc.'$type'.bw' ~/group/projects/vision/merged_input/200_noblack.11_22_2017.10kb.bed $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab'	
	time ~/group/software/ucsc/bigWigAverageOverBed $ct'.'$mk'.meanrc.'$type'.bw' ~/group/genome/mm10/mm10.1to19_X.genome.bed $ct'.'$mk'.meanrc.'$type'.sort.wg.tab'

	###
	sort -k1,1 $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab' > $ct'.'$mk'.meanrc.'$type'.sort.5kb.sort.tab' && mv $ct'.'$mk'.meanrc.'$type'.sort.5kb.sort.tab' $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab'
	sort -k1,1 $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab' > $ct'.'$mk'.meanrc.'$type'.sort.10kb.sort.tab' && mv $ct'.'$mk'.meanrc.'$type'.sort.10kb.sort.tab' $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab'

	### get bed files
	time cut -f1 ERY_ad.atac.meanrc.S3norm_NBP.sort.5kb.tab | awk -F '_' -v OFS='\t' '{print $1, $2, $3}' > $ct'.'$mk'.meanrc.bg.txt'
	### get 5kb & 10kb & wg
	time cut -f6 $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab' > $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab.txt'
	time cut -f6 $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab' > $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab.txt'
	time paste $ct'.'$mk'.meanrc.bg.txt' $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab.txt' $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab.txt' | sort -k1,1 -k2,2n > $ct'.'$mk'.meanrc.'$type'.sort.bg.tab.txt'
	### get MACS2 bg
	#time Rscript get_local_bg_rc.atac.R $ct'.'$mk'.meanrc.'$type'.sort.wg.tab' $ct'.'$mk'.meanrc.'$type'.sort.bg.tab.txt' $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph'
	time Rscript get_local_bg_rc.atac.R $ct'.'$mk'.meanrc.'$type'.txt' atac_commonpk.txt.cbg.txt $ct'.'$mk'.meanrc.'$type'.sort.bg.tab.txt' $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph'
	### rm tmp files
	#rm $ct'.'$mk'.meanrc.'$type'.sort.5kb.tab.txt' $ct'.'$mk'.meanrc.'$type'.sort.10kb.tab.txt'
	#rm $ct'.'$mk'.meanrc.bg.txt'
	#rm $ct'.'$mk'.meanrc.'$type'.sort.bg.tab.txt'	

	### MACS2 pk calling
#	time macs2 bdgcmp -t $ct'.'$mk'.meanrc.'$type'.sort.txt' -c $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' -m ppois -o $ct'.'$mk'.meanrc.'$type'.macs2pk.ppois.txt'
#	time macs2 bdgcmp -t $ct'.'$mk'.meanrc.'$type'.sort.txt' -c $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' -m qpois -o $ct'.'$mk'.meanrc.'$type'.macs2pk.qpois.txt'
#	time macs2 bdgcmp -t $ct'.'$mk'.meanrc.'$type'.sort.txt' -c $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' -m logLR -p 0.1 -o $ct'.'$mk'.meanrc.'$type'.macs2pk.logLR.txt'
	time Rscript get_NB_signal_track.R $ct'.'$mk'.meanrc.'$type'.sort.txt' $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' atac_commonpk.txt.cbg.sort.txt $ct'.'$mk'.meanrc.'$type'.sort'
	time Rscript get_NB_signal_track.R $ct'.'$mk'.meanrc.'$type'.sort.txt' atac_uniform_bg.txt atac_commonpk.txt.cbg.sort.txt $ct'.'$mk'.meanrc.'$type'.sort.ubg'
#	time macs2 bdgpeakcall -i $ct'.'$mk'.meanrc.'$type'.macs2pk.ppois.txt' -o $ct'.'$mk'.meanrc.'$type'.macs2pk.report.txt' --cutoff-analysis
done

