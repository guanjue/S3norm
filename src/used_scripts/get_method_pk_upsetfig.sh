declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

### get pooled peaks
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	ct=$(echo "$LINE" | awk '{print $3}')
	mk=$(echo "$LINE" | awk '{print $4}')
	echo $ct'.'$mk
	### prepare for method shared peaks
	cp $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.raw.macs2pk.NBQ.2.txt' \
	$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed'
	for type in "${type_vector[@]}"
	do
		### pool peaks
		cat $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed' \
		$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt' \
		> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed.tmp' \
		&& mv $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed.tmp' \
		$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed'
	done

	### get sort uniq pks
	time sort -k1,1 -k2,2n $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed' | uniq \
	> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.su.bed'
	time bedtools merge -i $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.su.bed' \
	> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.bed'
	time sort -k1,1 -k2,2n $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.bed' \
	> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed'

	### clean
	rm $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.bed'
	rm $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.bed'
	rm $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.su.bed'

	### intersect with method peaks
	for type in "${type_vector[@]}"
	do
		echo $type
		bedtools intersect -a $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed' \
		-b $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt' -c \
		> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt'
	done

	### cut columns with method peaks
	cp $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed' \
	$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt'
	for type in "${type_vector[@]}"
	do
		echo $type
		### cut
		cut -f4 $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt' \
		> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt.tmp'
		### paste
		paste $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt' \
		$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt.tmp' \
		> $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt.tmp' \
		&& mv $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt.tmp' \
		$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt'
		### clean data
		rm $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt.tmp'
		rm $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.'$type'.macs2pk.NBQ.2.txt.count.txt'
	done

	### plot upset figure
	for type in "${type_vector[@]}"
	do
		echo $type
		### plot
		time Rscript get_upsetfig.R $ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.pooled.merge.sort.bed.txt' \
		'upsetfig/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.upset.pdf'
	done
done < rc_list.atac.txt


