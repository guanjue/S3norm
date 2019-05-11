### get rnaHtseqCountsall matrix
cat rnaHtseqCountsall.merged.header.txt > rnaHtseqCountsall.merged.txt
tail -n+2 rnaHtseqCountsall.txt | awk -F '\t' -v OFS='\t' '{print $1, $2+$3, $6+$7, $8+$9, $10+$11, $12+$13, $14+$15, $16+$17, $18, $19+$20, $21+$22, $23+$24}' >> rnaHtseqCountsall.merged.txt

cut -f1,2,3,4 gencode_pc_sort.TSSexp2kb.bed > gencode_pc_sort.TSSexp2kb.coordi_id.bed
cut -f1,2,3,4 gencode_pc_sort.TSSexp5kb.bed > gencode_pc_sort.TSSexp5kb.coordi_id.bed
cut -f1,2,3,4 gencode_pc_sort.TSSexp10kb.bed > gencode_pc_sort.TSSexp10kb.coordi_id.bed

declare -a rna_ct_vector=("CFU_E_ad" "CMP" "ERY_ad" "GMP" "MK_imm_ad" "LSK_BM" "MEP" "MONO_BM" "NEU" "ER4" "G1E")
declare -a method_type=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")

### get pc h3k4me3 signal
for type in "${method_type[@]}"
do
	echo $type
	ct=${rna_ct_vector[1]}
	time ~/group/software/ucsc/bigWigAverageOverBed ../$ct'.h3k4me3.meanrc.'$type'.bw' gencode_pc_sort.TSSexp5kb.coordi_id.bed $ct'.h3k4me3.meanrc.'$type'.bw.tab'
	cut -f1 $ct'.h3k4me3.meanrc.'$type'.bw.tab' > 'h3k4me3.TSSexp5kb.'$type'.txt'
	### loop ct
	for ct in "${rna_ct_vector[@]}"
	do
		echo $type
		time ~/group/software/ucsc/bigWigAverageOverBed ../$ct'.h3k4me3.meanrc.'$type'.bw' gencode_pc_sort.TSSexp5kb.coordi_id.bed $ct'.h3k4me3.meanrc.'$type'.bw.tab'
		cut -f5 $ct'.h3k4me3.meanrc.'$type'.bw.tab' > $ct'.h3k4me3.meanrc.'$type'.bw.tab.txt'
		paste 'h3k4me3.TSSexp5kb.'$type'.txt' $ct'.h3k4me3.meanrc.'$type'.bw.tab.txt' > 'h3k4me3.TSSexp5kb.'$type'.txt.tmp' && mv 'h3k4me3.TSSexp5kb.'$type'.txt.tmp' 'h3k4me3.TSSexp5kb.'$type'.txt'
	done
	rm *'.h3k4me3.meanrc.'$type'.bw.tab.txt'
	rm *'.h3k4me3.meanrc.'$type'.bw.tab'
done

### select pc gene & sort rnaHtseqCountsall mat
time python vlookup.py -t rnaHtseqCountsall.merged.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o rnaHtseqCountsall.merged.pc.txt -k F
time python vlookup.py -t h3k4me3.TSSexp10kb.raw.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o h3k4me3.TSSexp10kb.raw.pc.txt -k F
time python vlookup.py -t h3k4me3.TSSexp10kb.TSnorm.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o h3k4me3.TSSexp10kb.TSnorm.pc.txt -k F
time python vlookup.py -t h3k4me3.TSSexp10kb.MAnorm.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o h3k4me3.TSSexp10kb.MAnorm.pc.txt -k F
time python vlookup.py -t h3k4me3.TSSexp10kb.QTnorm.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o h3k4me3.TSSexp10kb.QTnorm.pc.txt -k F
time python vlookup.py -t h3k4me3.TSSexp10kb.S3norm.txt -m 1 -s gencode.vM4.annotation.pc.sorted.bed -n 4 -o h3k4me3.TSSexp10kb.S3norm.pc.txt -k F






