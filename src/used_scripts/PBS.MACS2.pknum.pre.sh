#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/scratch/vision/all_final_data/rc_norm

declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "POIS" "Z" "S3norm" "NBP" "S3norm_NBP")

while read LINE
do
sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
mkdir $ct'_'$mk'_macs2pk'
### get pk calling result and pk number
for type in "${type_vector[@]}"
do
	rm $ct'.'$mk'.pknum_vs_thresh.'$type'.txt'
	echo $type
	for t in $(seq 0 5 100)
	do
		echo $t
		time macs2 bdgpeakcall -i $ct'.'$mk'.meanrc.'$type'.sort.txt' -o $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' -c $t
		wc -l $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' >> $ct'.'$mk'.pknum_vs_thresh.'$type'.txt'
		mv $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' $ct'_'$mk'_macs2pk/'	
	done
done
done < rc_list.atac.txt

