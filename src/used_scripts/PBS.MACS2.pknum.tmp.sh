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
file=$(head -${PBS_ARRAYID} rc_list.atac.txt | tail -1)
echo $file

#declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "POIS" "Z" "S3norm" "NBP" "S3norm_NBP")
declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")
declare -a type_vector=("S3norm")

LINE=$(head -${PBS_ARRAYID} rc_list.atac.txt | tail -1)
sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $sig1 
echo $sig2
echo $ct
echo $mk
mkdir $ct'_'$mk'_macs_pk'
### for RC
for type in "${type_vector[@]}"
do
	echo $type
	for t in $(seq 0 1 100)
	do
		echo $t
		### POIS p
#		time macs2 bdgcmp -t $ct'.'$mk'.meanrc.'$type'.sort.txt' -c $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' -m ppois -o $ct'.'$mk'.meanrc.'$type'.macs2pk.ppois.txt'
#		time macs2 bdgpeakcall -i $ct'.'$mk'.meanrc.'$type'.macs2pk.ppois.txt' -o $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' -c $t
#		mv $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' $ct'_'$mk'_macs_pk/'
		### POIS q
		time macs2 bdgcmp -t $ct'.'$mk'.meanrc.'$type'.sort.txt' -c $ct'.'$mk'.meanrc.'$type'.sort.bg.bedgraph' -m qpois -o $ct'.'$mk'.meanrc.'$type'.macs2pk.qpois.txt'
		time macs2 bdgpeakcall -i $ct'.'$mk'.meanrc.'$type'.macs2pk.qpois.txt' -o $ct'.'$mk'.meanrc.'$type'.macs2pk.qpois.'$t'.txt' -c $t
		mv $ct'.'$mk'.meanrc.'$type'.macs2pk.qpois.'$t'.txt' $ct'_'$mk'_macs_pk/'
	done
done

### for POIS and NBP
declare -a type_vector_p=("POIS" "NBP" "S3norm_NBP")
for type in "${type_vector_p[@]}"
do
	echo $type
	for t in $(seq 0 1 100)
	do
		echo $t
		time macs2 bdgpeakcall -i $ct'.'$mk'.meanrc.'$type'.sort.txt' -o $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' -c $t
		mv $ct'.'$mk'.meanrc.'$type'.macs2pk.'$t'.txt' $ct'_'$mk'_macs_pk/'
	done
done

