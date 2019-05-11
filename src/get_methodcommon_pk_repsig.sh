declare -a method_type=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

for type in "${method_type[@]}"
do
echo $type
while read LINE
do
sig1=$(echo "$LINE" | awk '{print $1}')
sig2=$(echo "$LINE" | awk '{print $2}')
ct=$(echo "$LINE" | awk '{print $3}')
mk=$(echo "$LINE" | awk '{print $4}')
echo $ct'.'$mk
### get signal for method shared peaks
#wc -l ../rc_norm/$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed'
### get signal
time ~/group/software/ucsc/bigWigAverageOverBed $ct'rep1.'$mk'.meanrc.'$type'.bw' ../rc_norm/$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' $ct'rep1.'$mk'.meanrc.'$type'.bw.tab'
time ~/group/software/ucsc/bigWigAverageOverBed $ct'rep2.'$mk'.meanrc.'$type'.bw' ../rc_norm/$ct'_'$mk'_macs_pk/'$ct'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' $ct'rep2.'$mk'.meanrc.'$type'.bw.tab'
### cut columns
cut -f5 $ct'rep1.'$mk'.meanrc.'$type'.bw.tab' > $ct'rep1.'$mk'.meanrc.'$type'.bw.tab.txt'
cut -f5 $ct'rep2.'$mk'.meanrc.'$type'.bw.tab' > $ct'rep2.'$mk'.meanrc.'$type'.bw.tab.txt'
paste $ct'rep1.'$mk'.meanrc.'$type'.bw.tab.txt' $ct'rep2.'$mk'.meanrc.'$type'.bw.tab.txt' > $ct'rep12.'$mk'.meanrc.'$type'.bw.tab.txt'
### rm 
rm $ct'rep1.'$mk'.meanrc.'$type'.bw.tab.txt'
rm $ct'rep2.'$mk'.meanrc.'$type'.bw.tab.txt'
rm $ct'rep1.'$mk'.meanrc.'$type'.bw.tab'
rm $ct'rep2.'$mk'.meanrc.'$type'.bw.tab'
done < rc_list.atac.norep.txt
done
