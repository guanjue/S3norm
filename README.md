# S3norm
## From  Reads Count (RC) to S3norm normalized -log10(p-value)

#### Motivation: The quantitative comparison of epigenomic data across multiple cell types has become a promising way to understand the biological function of epigenetic modifications. Due to difference in sequencing depth and signal-to-noise ratio, however, the raw epigenomic data may not reflect the real biological difference between cell types. Existing normalization methods are mainly designed for scaling signals in either the whole-genome or the peak regions, without considering the potentially different scaling factors between peak and background regions. Results: We propose a new data normalization method, S3norm, that normalizes the data by using a monotonic nonlinear data transformation to match signals in both the peak regions and the background regions differently, such that both sequencing depth and signal-to-noise ratio between data sets can be simulatenously normalized. We show that the S3norm normalized epigenomic data can better reflect real biological differences across multiple cell types.



<img src="https://github.com/guanjue/S3norm/blob/master/example_figures/overall_pipeline.png" width="800"/>

##### Figure 1. The overall workflow of the S3norm normalization method. There are three major steps in S3norm. (a) The 1st step is to convert reads count in the 200-bp bins to -log10(p-value) for each epigenomic dataset. Each box represents the different signal tracks of the same data. The first one is the raw reads count of the G1E H3K4me3 dataset. The second one is the reads count of input sample. The third one is the -log10(p-value) of the G1E H3K4me3 dataset. The shoulder of the peaks are reduced after convert reads count to -log10(p-value) with the background adjustment. (b) The 2nd step is selecting the dataset with the highest SNR as the reference dataset for the S3norm normalization. The barplot represents the SNRs of all datasets. The dataset with the highest SNR (dataset with the orange bar) will be selected as the reference dataset. (c) The 3rd step is using a monotonic nonlinear data transformation model to normalize both the SNR and SD between the two datasets. The (1) part is identifing common peak regions and the common background regions between the two datasets. The left scatterplot is showing the signal of each bin in the target dataset and the reference dataset. In the right scatterplot, each data point is colored based on the type of the data point. The orange data points represent the common peak bins. The gray data points represent the common background bins. The blue data points represent the dataset-specific bins. The (2) part is using the monotonic nonlinear data transformation model to rotate the signal of the target dataset, so that (i) the means of the common peak regions of two datasets and (ii) the means of the common background regions of the two datasets can be matched. The details of each step is discussed in the method section.




## Install S3norm pipeline
#### Clone the github repository 
```
git clone https://github.com/guanjue/S3norm.git
```
#### Install dependency: go to S3norm folder and run the following command
```
time bash INSTALL.sh
```



#### The input file list for S3norm: 
##### The file name should contain the cell type name and the mark name and sample id separated by ".":
###### cell_type.mark_name.sample_id.canb_be_anything
###### e.g. B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal

##### If the user want to keep the replicates separate, the 'cell_type' in the file name should be replaced by the 'cell_type_sample_id'
###### cell_type_sample_id.mark_name.canb_be_anything
###### e.g. B_SPL_100035.h3k27acrep.bamtobed5endintersect.signal

###### 1st column: the signal of each bin;
###### all of the input file should be saved in the input folder (input_dir)
```
>>> head -1000000 B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal | tail -20
0
0
0
1
1
0
0
0
2
0
0
0
0
4
```

##### The input filename list for S3norm: each column is separated by tab
###### 1st column: target dataset; 
###### 2nd column: the no antibody control file for the target dataset; Here, we used the same merge input file base on the 21 input files from 11 cell types. For each input file, it is first normalized to ER4.137.input.signal input dataset based the ratio of total reads count. Then, the same merge input file is the mean signal of the 21 normalized the input files. 
###### For the atac-seq without input signal. We used a input file with all bins equal to one as the input signal.
```
info_table_all.rc2nbp.txt
>>> head info_table_all.rc2nbp.txt
B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CFU_E_ad.h3k27acrep.100030.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CLP.h3k27acrep.100039.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CMP.h3k27acrep.100027.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.538.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.539.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
......
```


## Run S3norm
##### (1) copy the 'run_pipeline.sh' in the S3norm into the working directory
```
cp ~/group/software/S3norm/run_pipeline.sh working_dir/
```
##### (2) change the following parameters in the 'run_IDEAS.sh' file:
###### script_dir='absolute path to the IDEAS_2018 dir'
###### working_dir='absolute path to the working directory'
###### input_file_list='the input and the correponding no antibody control file list'
###### input_dir='the input file's folder'
###### overall_upper='the upper limit of the output file'
###### overall_lower='the lower limit of the output file'
###### select_method='method used to select reference dataset (frip/snr)'
###### select_ref_version='method used to select reference dataset (max/median)'
###### user_given_global_ref='user given global reference dataset (if empty, pipeline will user the dataset with the highest frip/snr score dataset)'
###### bin_num='number of bins'
```
>>> head -100 run_IDEAS.sh 
###### set parameters
script_dir=/storage/home/gzx103/group/software/S3norm/src/
working_dir=/storage/home/gzx103/scratch/S3norm/test_pipeline/
input_dir=/storage/home/gzx103/scratch/S3norm/test_pipeline/input_5end_rc/
input_file_list=info_table_all.rc2nbp.txt
overall_upper=100
overall_lower=0
select_method=frip
select_ref_version=median
user_given_global_ref=NA
bin_num=13554672

time bash $script_dir'overall_pipeline.sh' $script_dir $working_dir $input_dir $input_file_list $overall_upper $overall_lower $select_method $select_ref_version $user_given_global_ref $bin_num
```
##### (3) create the input_list_file
```
>>> head info_table_all.rc2nbp.txt
B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CFU_E_ad.h3k27acrep.100030.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CLP.h3k27acrep.100039.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CMP.h3k27acrep.100027.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.538.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.539.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
```

##### (4) use 'run_IDEAS.sh' script to run S3norm pipeline
```
time bash run_pipeline.sh
```



## Output results for test data
### All output files will be saved to the following directories inside the working directory:
```
nbp/: The p-value of each sample based on the NB background model.
s3norm_info/: The scatterplot and the parameters used in S3norm normalization for each dataset
s3norm_sig/: The signal files of the S3norm normalized data
ref_info/: The scatterplot and the parameters used in S3norm normalization for reference datasets across all marks
s3norm_ref_sig/: The signal files of the S3norm normalized reference data
fisherp/: The merged p-value from each sample's NB p-value by the Fisher's method
list_files/: All of the list files used in the pipeline
s3norm_0_100_sig/: The signal files of the S3norm normalized data with the upper & lower bound limitation
```


## References

##### S3norm & vision paper

