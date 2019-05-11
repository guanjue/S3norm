declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	ct=$(echo "$LINE" | awk '{print $3}')
	mk=$(echo "$LINE" | awk '{print $4}')
	echo $ct'.'$mk
	### prepare for method shared peaks
	cp $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.raw.macs2pk.NBQ.2.txt' \
	$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed'
	for type in "${type_vector[@]}"
	do
		### intersect with peaks
		bedtools intersect -a $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' \
		-b $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt' \
		-wa -u > $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed.tmp'
		### change filename
		mv $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed.tmp' \
		$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed'
	done
	### check shared peak number
	wc -l $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed'
done < rc_list.atac.txt

