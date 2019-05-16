mkdir sigdif_pk

declare -a type_vector=("raw" "TSnorm" "MAnorm" "QTnorm" "S3norm")
declare -a rep_vector=("rep1" "rep2")
declare -a ct1_vector=("ERY_ad" "CFU_E_ad" "CFUMK" "CMP" "ER4" "G1E" "GMP" "LSK_BM" "MEP" "NEU")
declare -a ct2_vector=("ERY_ad" "CFU_E_ad" "CFUMK" "CMP" "ER4" "G1E" "GMP" "LSK_BM" "MEP" "NEU")

#ct1=CMP
#ct2=G1E
mk=atac

for ct1 in "${ct1_vector[@]}"
do
	for ct2 in "${ct2_vector[@]}"
	do
		### merge two cell type peaks
		cat $ct1'_'$mk'_macs_pk/'$ct1'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' \
		$ct2'_'$mk'_macs_pk/'$ct2'.'$mk'.meanrc.macs2pk.NBQ.2.shared.bed' \
		| sort -k1,1 -k2,2n > 'sigdif_pk/'$ct1'.'$ct2'.shared.bed'
		### merge peak
		bedtools merge -i 'sigdif_pk/'$ct1'.'$ct2'.shared.bed' > 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.bed'
		cat 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.bed' | \
		awk -F '\t' -v OFS='\t' '{print $1,$2,$3, $1"_"$2"_"$3}' > 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed'
		###### get signal
		### type
		# ct1 sig
		for type in "${type_vector[@]}"
		do
			### get pkid
			time ~/group/software/ucsc/bigWigAverageOverBed \
			../rc_norm_rep/$ct1'rep1.'$mk'.meanrc.'$type'.bw' \
			'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed' \
			'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.rep1.'$type'.tab'
			time cut -f1 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.rep1.'$type'.tab' \
			> 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt'
			for ri in "${rep_vector[@]}"
			do
				time ~/group/software/ucsc/bigWigAverageOverBed \
				../rc_norm_rep/$ct1$ri'.'$mk'.meanrc.'$type'.bw' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.'$ri'.'$type'.tab'
				# ct2 sig
				time ~/group/software/ucsc/bigWigAverageOverBed \
				../rc_norm_rep/$ct2$ri'.'$mk'.meanrc.'$type'.bw' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct2'.'$ri'.'$type'.tab'
				### cut sig ct1
				time cut -f5 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.'$ri'.'$type'.tab' \
				> 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.'$ri'.'$type'.tab.txt'
				### paste to matrix
				time paste 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct1'.'$ri'.'$type'.tab.txt' \
				> 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt.tmp' && \
				mv 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt.tmp' 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt'
				### cut sig ct2
				time cut -f5 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct2'.'$ri'.'$type'.tab' \
				> 'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct2'.'$ri'.'$type'.tab.txt'
				### paste to matrix
				time paste 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt' \
				'sigdif_pk/'$ct1'.'$ct2'.shared.merged.withid.bed.'$ct2'.'$ri'.'$type'.tab.txt' \
				> 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt.tmp' && \
				mv 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt.tmp' 'sigdif_pk/'$ct1'.'$ct2'.repsig.'$type'.txt'
			done
		done
		### mv file to folder
		mkdir sigdif_pk/$ct1'_'$ct2'_sharedpk/'
		mv 'sigdif_pk/'$ct1'.'$ct2* sigdif_pk/$ct1'_'$ct2'_sharedpk/'
	done
done


