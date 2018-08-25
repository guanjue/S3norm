###### read parameters from inputs
script_dir=$1
working_dir=$2
input_dir=$3
input_file_list=$4
overall_upper=$5
overall_lower=$6
select_method=$7
select_ref_version=$8
user_given_global_ref=$9
bin_num=$10


### set rank lim & plot point number
rank_lim=$((bin_num / 100))
plot_num=$((bin_num / 20))

echo $rank_lim
echo $plot_num

### select top reference dataset for cross mark s3norm
if time Rscript $script_dir'get_top_ref.R' '.ref_frip.txt' $select_method $working_dir cross_mark_ref_list.txt $select_ref_version $user_given_global_ref; then echo 'select top reference dataset for cross mark s3norm DONE'; else echo 'ERROR: select top reference dataset for cross mark s3norm' && exit 1; fi



###### s3norm normalize reference datasets of all marks
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
	upperlim=323
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
	### peak norm
	if time python $script_dir's3norm.py' -r $sig1'.upperlim.txt' -t $sig2'.upperlim.txt' -m 1 -i 2 -f 0.05 -n $plot_num -l $rank_lim -a 100 -b 0 -s $script_dir -p p -c T; then echo 's3norm across datasets DONE'; else echo 'ERROR: s3norm across datasets' && exit 1; fi
	### rm tmp files
	rm $sig1'.upperlim.txt'
	rm $sig2'.upperlim.txt'
done < cross_mark_ref_list.txt.info.txt

### move ref norm files into ref_info folder
if [ -d $working_dir'ref_info/' ]; then echo $working_dir'ref_info/' exist; else mkdir $working_dir'ref_info/'; fi
mv *.s3norm.scatterplot.png $working_dir'ref_info/'
mv *.scatterplot.png $working_dir'ref_info/'
mv *.ref.info.txt $working_dir'ref_info/'



###### s3norm across datasets with the same mark
for mk in $(cat mark_list.txt)
do
	echo $mk
	while read LINE
	do
		sig1=$(echo "$LINE" | awk '{print $1}')
		sig2=$(echo "$LINE" | awk '{print $2}')
		sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
		upperlim=323
		lowerlim=0
		echo $sig1 
		echo $sig2
		echo $sig2_celltype
		### set upper limit
		cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
		cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
		### peak norm
		if time python $script_dir's3norm.py' -r $sig1'.upperlim.txt' -t $sig2'.upperlim.txt' -m 1 -i 2 -f 0.05 -n $plot_num -l $rank_lim -a 100 -b 0 -s $script_dir -p p -c F; then echo 's3norm across datasets DONE'; else echo 'ERROR: s3norm across datasets' && exit 1; fi
		### rm tmp files
		rm $sig1'.upperlim.txt'
		rm $sig2'.upperlim.txt'
	done < $mk'.ref_frip.txt.info.txt'
done
### move s3norm files into s3norm_sig folder
if [ -d $working_dir's3norm_info/' ]; then echo $working_dir's3norm_info/' exist; else mkdir $working_dir's3norm_info/'; fi
mv *.s3norm.scatterplot.png $working_dir's3norm_info/'
mv *.scatterplot.png $working_dir's3norm_info/'
mv *.info.txt $working_dir's3norm_info/'



###### mv S3norm normalized signal files & unnormalized signal files into *_sig folders
### s3norm signal
if [ -d $working_dir's3norm_sig/' ]; then echo $working_dir's3norm_sig/' exist; else mkdir $working_dir's3norm_sig/'; fi
mv *.s3norm.txt $working_dir's3norm_sig/'
### ref s3norm signal
if [ -d $working_dir's3norm_ref_sig/' ]; then echo $working_dir's3norm_ref_sig/' exist; else mkdir $working_dir's3norm_ref_sig/'; fi
mv *.s3norm.ref.txt $working_dir's3norm_ref_sig/'
### ref frip & snp
mv *.ref_frip.txt $working_dir'ref_info/'
### fisher pvalue signal without normalization
if [ -d $working_dir'fisherp/' ]; then echo $working_dir'fisherp/' exist; else mkdir $working_dir'fisherp/'; fi
mv *.fisher_p.txt $working_dir'fisherp/'
mv *.frip_snr.txt $working_dir'fisherp/'



###### set limit for signals
for filename in $(cat cell_marker_list.txt)
do
	echo $filename
	if cat $working_dir's3norm_sig/'$filename'.s3norm.txt' | awk -F '\t' -v OFS='\t' -v ul=$overall_upper -v ll=$overall_lower '{if ($1<ll) print ll; else if ($1>ul) print ul; else print $1}' > $filename'.s3norm.'$overall_lower'_'$overall_upper'.txt'; then echo 'set limit for signals DONE'; else echo 'ERROR: set limit for signals' && exit 1; fi
done
### mv to the output folder
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'



### list files
if [ -d $working_dir'list_files/' ]; then echo $working_dir'list_files/' exist; else mkdir $working_dir'list_files/'; fi
mv *list.txt $working_dir'list_files/'





