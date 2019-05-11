declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	ct=$(echo "$LINE" | awk '{print $3}')
	mk=$(echo "$LINE" | awk '{print $4}')
	echo $ct'.'$mk
	#wc -l $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed'
	cat $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.pkid.bed'
done < rc_list.atac.txt


