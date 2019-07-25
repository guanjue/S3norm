# S3norm: 
## Simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data
### From  Reads Count (RC) to S3norm normalized -log10(p-value)

#### Quantitative comparison of epigenomic data across multiple cell types or experimental conditions is a promising way to understand the biological functions of epigenetic modifications. Difference in sequencing depth and signal-to-noise ratios resulted from different experiments however hinder our ability to identify real biological variation from raw epigenomic data. Proper normalization is therefore required prior to data analyses to gain meaningful insights. Existing data normalization methods mostly standardize signals by rescaling either background regions or peak regions, assuming that the same scale factor is applicable to both background regions and peak regions. While such methods adjust for differences due to sequencing depths, they fail when the signal-to-noise ratios are different across experiments. We propose a new data normalization method, called S3norm, that normalizes the sequencing depths and signal-to-noise ratios across different data sets simultaneously by a monotonic nonlinear transformation. Empirically we show that the epigenomic data normalized by our method can better capture real biological variation, such as their impact one gene expression regulation and the numbers of the their occurrences, across different cell types than existing methods. 



<img src="https://github.com/guanjue/S3norm/blob/master/example_figures/overall_pipeline.png" width="800"/>

##### Figure 1. Overview of the S3norm method. These are the scatterplots of read counts (log scale) in 10,000 randomly selected genome locations (200bp) in target cell (x-axis) and reference cell (y-axis). The left figure is the signal before S3norm. The right figure is the signal after S3norm. The S3norm applies a monotonic nonlinear model ( log⁡(Y_(norm,i))=log(α)+βlog(Y_i) ) to rotate the target signal so that (1) the mean signals of common peaks (green point, highlighted by black dash circle) and (2) the mean signals of common background (dark blue point, highlighted by black dash circle)  can be matched between the two data sets. The original data were split into three groups: the common peak regions (orange), the common background regions (gray), and the rest bins (blue). The overall mean is represented by a black point. 



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
####
#### !!! The first three columns of All of bedgraph files in the filelist should be exactly the same !!!. 
#### !!! Only the fourth column is different !!!. 
####
##### Each of bedgraph file in the filelist (both ChIP bedgraph and Control bedgraph) need to be sorted before running S3norm. 
##### This can be done by the following command:
```
>>> sort -k1,1 -k2,2n sig1.UNsorted.bedgraph > sig1.sorted.bedgraph
```


## Run S3norm
### (1) User can run the full pipeline of S3norm by using 'S3norm_pipeline.py' scripts in the 'S3norm/src/' folder
##### Before running S3norm full pipeline, user should let S3norm know where is the work directory and where is the script directory.
##### This can be done by replacing the "script_directory='/Users_given_directory/S3norm/'" 
##### and the "working_script_directory='/Users/universe/software/S3norm/example_file/'" 
##### !!! All of the bedgraph file should in the working directory !!!
```
### Setting script directory
script_directory='/Users/universe/software/S3norm/'

### Setting working directory
working_script_directory='/Users/universe/software/S3norm/example_file/'

### Entering working directory
cd $working_script_directory

##################
### Run S3norm
##################
time python $script_directory'/src/S3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt

```

### (2) User can also run each step in S3norm by using the following scripts in the 'S3norm/src/' folder





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

