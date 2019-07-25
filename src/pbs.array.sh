module load gcc/5.3.1 
module load python/2.7.14-anaconda5.0.1
module load bedtools
module load r

### after get_NBP.sh
script_dir='/storage/work/bmg137/ideas/S3norm/src/'
working_dir='/gpfs/scratch/bmg137/ideasBp19/'

ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1"."$2}' | sort -u > cell_marker_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $2}' | sort -u > mark_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1}' | sort -u > cell_list.txt
mv *.nbp_2r_bgadj.txt $working_dir'nbp/'
mv *.mvsp.txt $working_dir'nbp/'


### after get_fisherp.sh

### after get_median_ref.sh

### after get_ref_bw.sh

script_dir='/storage/work/bmg137/ideas/S3norm/src/'
working_dir='/gpfs/scratch/bmg137/ideasBp20/'

select_method=frip


for mk in $(cat mark_list.txt)
do
	if [ -z "$mk" ]; then continue; fi
	echo "Doing mark $mk" 1>&2
	ls *$mk*.frip_snr.txt > $mk'.file_list.txt'
	echo $select_method 1>&2
	time Rscript 'get_mk_ref.R' $mk'.file_list.txt' $select_method $mk'.ref_'$select_method'.txt' $mk'.fisherp.ref.txt' $mk'.fisherp.ref.info.txt' $mk'.fisherp.s3norm.ref.txt'
#	if time Rscript 'get_mk_ref.R' $mk'.file_list.txt' frip $mk'.ref_'$select_method'.txt'; then echo 'select reference dataset for s3norm'; else echo 'ERROR: select reference dataset for s3norm' && exit 1; fi
done

time Rscript $script_dir'get_top_ref.R' '.ref_'$select_method'.txt' snr $working_dir cross_mark_ref_list.txt max NA


LINE=$(head -1 cross_mark_ref_list.txt.info.txt | tail -1)
sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
upperlim=100
lowerlim=0
echo $sig1 
echo $sig2
echo $sig2_celltype
time cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
time cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt'
time python $script_dir's3norm.py' -r $sig1'.upperlim.txt' -t $sig2'.upperlim.txt' -m 1 -i 2 -f 0.05 -n 10000 -l 10000 -a 100 -b 0 -s $script_dir -p z -c T

### after get_s3normref.sh
for mk in $(cat mark_list.txt)
do
	echo $mk
	head -15 $mk'.ref_'$select_method'.txt.info.txt' > $mk'.ref_'$select_method'.txt.info.1.txt'
	tail -n+16 $mk'.ref_'$select_method'.txt.info.txt' > $mk'.ref_'$select_method'.txt.info.2.txt'
done


script_dir='/storage/work/bmg137/ideas/S3norm/src/'
working_dir='/gpfs/scratch/bmg137/ideasBp19/'
	
if [ -d $working_dir'ref_info/' ]; then echo $working_dir'ref_info/' exist; else mkdir $working_dir'ref_info/'; fi
#mv *.s3norm.scatterplot.png $working_dir'ref_info/'
#mv *.scatterplot.png $working_dir'ref_info/'
mv *.ref.info.txt $working_dir'ref_info/'

### after get_upperlimed.sh

### after get_upperlimedref.sh

### get cbg cbp get_cpkbg_p.sh

### after get_s3norm.1.sh
### after get_s3norm.2.sh
script_dir='/storage/work/bmg137/ideas/S3norm/src/'
working_dir='/gpfs/scratch/bmg137/ideasBp20/'

if [ -d $working_dir's3norm_info/' ]; then echo $working_dir's3norm_info/' exist; else mkdir $working_dir's3norm_info/'; fi
#mv *.s3norm.scatterplot.png $working_dir's3norm_info/'
#mv *.scatterplot.png $working_dir's3norm_info/'
mv *.info.txt $working_dir's3norm_info/'


###### mv S3norm normalized signal files & unnormalized signal files into *_sig folders
### s3norm signal
if [ -d $working_dir's3norm_sig/' ]; then echo $working_dir's3norm_sig/' exist; else mkdir $working_dir's3norm_sig/'; fi
mv *.s3norm.txt $working_dir's3norm_sig/'
### ref s3norm signal
if [ -d $working_dir's3norm_ref_sig/' ]; then echo $working_dir's3norm_ref_sig/' exist; else mkdir $working_dir's3norm_ref_sig/'; fi
mv *.s3norm.ref.txt $working_dir's3norm_ref_sig/'
### ref frip & snp
mv *'.ref_'$select_method'.txt' $working_dir'ref_info/'
### fisher pvalue signal without normalization
if [ -d $working_dir'fisherp/' ]; then echo $working_dir'fisherp/' exist; else mkdir $working_dir'fisherp/'; fi
mv *.fisher_p.txt $working_dir'fisherp/'
mv *.frip_snr.txt $working_dir'fisherp/'

### after get_upperlim_lowerlim_output.sh
### mv to the output folder
overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'

### after get_upperlim_lowerlim_output_2_100.sh
overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig/'

### after get_upperlim_lowerlim_output_2_100.sh
overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_max/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_max/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_max/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_max/'

### after get_upperlim_lowerlim_output_2_100.sh
overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median/'

### after get_upperlim_lowerlim_output_2_100.sh
overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_2/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_3/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_4/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_5/'


overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_6/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_7/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_8/'


overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'

overall_lower=3
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_9/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_11/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_12/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_13/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_14/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_15/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_16/'


overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_17/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_18/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_19/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_20/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_21/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_22/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23a/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23a/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23a/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23a/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_23/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_24/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.a.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.a.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_25a/'

overall_lower=0
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_26/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1/'

overall_lower=1
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.a.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/'

overall_lower=2
overall_upper=100
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.a.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1a/'

overall_lower=0
overall_upper=2
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.b.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/'

overall_lower=0
overall_upper=1
if [ -d $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/' ]; then echo $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/' exist; else mkdir $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/'; fi
mv *'.s3norm.'$overall_lower'_'$overall_upper'.b.txt' $working_dir's3norm_'$overall_lower'_'$overall_upper'_sig_median_1b/'

### list files
if [ -d $working_dir'list_files/' ]; then echo $working_dir'list_files/' exist; else mkdir $working_dir'list_files/'; fi
mv *list.txt $working_dir'list_files/'





