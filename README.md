# S3norm
## From  Reads Count (RC) to S3norm normalized -log10(p-value)

#### Motivation: The quantitative comparison of epigenomic data across multiple cell types has become a promising way to understand the biological function of epigenetic modifications. Due to difference in sequencing depth and signal-to-noise ratio, however, the raw epigenomic data may not reflect the real biological difference between cell types. Existing normalization methods are mainly designed for scaling signals in either the whole-genome or the peak regions, without considering the potentially different scaling factors between peak and background regions. Results: We propose a new data normalization method, S3norm, that normalizes the data by using a monotonic nonlinear data transformation to match signals in both the peak regions and the background regions differently, such that both sequencing depth and signal-to-noise ratio between data sets can be simulatenously normalized. We show that the S3norm normalized epigenomic data can better reflect real biological differences across multiple cell types.



<img src="https://github.com/guanjue/S3norm/blob/master/example_figures/overall_pipeline.png" width="800"/>

##### Figure 1. The overall workflow of the S3norm normalization method. There are three major steps in S3norm. (a) The 1st step is to convert reads count in the 200-bp bins to -log10(p-value) for each epigenomic dataset. Each box represents the different signal tracks of the same data. The first one is the raw reads count of the G1E H3K4me3 dataset. The second one is the reads count of input sample. The third one is the -log10(p-value) of the G1E H3K4me3 dataset. The shoulder of the peaks are reduced after convert reads count to -log10(p-value) with the background adjustment. (b) The 2nd step is selecting the dataset with the highest SNR as the reference dataset for the S3norm normalization. The barplot represents the SNRs of all datasets. The dataset with the highest SNR (dataset with the orange bar) will be selected as the reference dataset. (c) The 3rd step is using a monotonic nonlinear data transformation model to normalize both the SNR and SD between the two datasets. The (1) part is identifing common peak regions and the common background regions between the two datasets. The left scatterplot is showing the signal of each bin in the target dataset and the reference dataset. In the right scatterplot, each data point is colored based on the type of the data point. The orange data points represent the common peak bins. The gray data points represent the common background bins. The blue data points represent the dataset-specific bins. The (2) part is using the monotonic nonlinear data transformation model to rotate the signal of the target dataset, so that (i) the means of the common peak regions of two datasets and (ii) the means of the common background regions of the two datasets can be matched. 



## Require python 2.7 and R

## Install S3norm pipeline
#### Clone the github repository 
```
git clone https://github.com/guanjue/S3norm.git
```
#### Install dependency: go to S3norm folder and run the following command
```
time bash INSTALL.sh
```

## Input files for S3norm
### S3norm uses bedgraph files as input files
#### (1) The input filelist for S3norm
##### The filelist contains the names of the ChIP bedgraph file and the control bedgraph. Each row is one ChIP-seq sample. The 1st column is the ChIP bedgraph and 2nd column is the Contrl bedgraph file. (Separated by tab "\t") 
##### The example of the filelist is in the 'example_file' folder.
```
>>> head file_list.txt 
sig1.bedgraph	sig1.ctrl.bedgraph
sig2.bedgraph	sig2.ctrl.bedgraph
sig3.bedgraph	sig3.ctrl.bedgraph
```

#### (2) The bedgraph files for S3norm
##### The format of bedgraph is as follows:
##### The four columns are chromosome, bin_start, bin_end, and signal in the bin.
##### For the S3norm full pipeline, the average read counts for each bin should be used as signal.
```
>>> head sig1.bedgraph
chrX	23515400	23515600	30.95
chr10	97283000	97283200	83.61
chr9	82643200	82643400	56.96
chr4	1898800	1899000	115.65
chrX	155824800	155825000	65.42
chr8	100720200	100720400	444.49
chrY	2919000	2919200	0
chr6	25794000	25794200	97.7
chr4	190032400	190032600	111.24
chr17	7828000	7828200	76.11
```

#### !!! The first three columns of All of bedgraph files in the filelist should be exactly the same !!!. 
#### !!! Only the fourth column is different !!!. 
####
##### Each of bedgraph file in the filelist (both ChIP bedgraph and Control bedgraph) need to be sorted before running S3norm. 
##### This can be done by the following command:
```
>>> sort -k1,1 -k2,2n sig1.UNsorted.bedgraph > sig1.sorted.bedgraph
```


## Run S3norm
### (1) Use 'S3norm_pipeline.py' to run S3norm pipeline
```
### Setting script directory
script_directory='/Users/universe/Documents/2018_BG/S3norm/'
### Setting working directory
working_script_directory='/Users/universe/Documents/2018_BG/S3norm/example_file'
### Entering working directory
cd $working_script_directory
### Run S3norm
time python $script_directory'/src/S3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt
```


### Options
#### Required parameters:
##### For S3norm, there are just two required parameters.
##### (1) -s script_folder: the script directory of S3norm (e.g. /Users/universe/Documents/2018_BG/S3norm/src/)
##### (2) -s input_file_list: the input filelist for S3norm (e.g. The 'file_list.txt' in the 'S3norm/example_file/' folder)

#### For other parameters, user can check them by using the following command:
```
python $script_directory'/src/S3norm_pipeline.py' -s script_folder -t input_file_list -r (reference_method: mean or median) -m (Method for matching peaks and background: non0mean, non0median, mean, median) -i initial_B -f FDR_thresh -l rank_lim_p -a upperlimit -b lowerlimit -p (p-value_method: neglog10p, z) -k common_pk_binary (0 for nocommon_pk; common_pk_binary.txt) -g common_bg_binary (0 for nocommon_pk; common_bg_binary.txt)
```

#### The other parameters that can be changed:
```
##### (1) -r : The method for generating the reference signal track. Options: median (default) or mean
##### (2) -m : The method for matching peaks and background. Options: non0mean (default), non0median, mean, median)
##### (3) -i : The initial value for the power parameter in the non-linear transformation. Default: 2.0
##### (4) -f : The FDR threshold for identifying common peaks. Default: 0.05 . The range is 0.0 < x < 1.0
##### (5) -l : The minimum proportion of bins are used as peak for S3norm. Default: 0.001 . The range is 0.0 < x < 1.0
##### (6) -a : The upperlimit for signal. This is to reduce the bias cause by extrame signals in the data. Default: 100000
##### (7) -b : The lowerlimit for signal. S3norm requires all signal to be x >= 0 . Default: 0
##### (8) -p : The method used to identify common peaks. Options: z (Default) and neglog10p (negative log10 p-value from background model)
##### (9) -k : The user given common peak regions. Options: 0 (Default, the common peak will be identified by S3norm) or filename (a file points out which bins are the common peaks. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common peak, the column should be 1 for that row. Otherwise, it should be 0 )
##### (10) -g  : The user given common background regions. Options: 0 (Default, the common background will be identified by S3norm) or filename (a file points out which bins are the common background. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common background, the column should be 1 for that row. Otherwise, it should be 0 )
```









## Output results for S3norm
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

