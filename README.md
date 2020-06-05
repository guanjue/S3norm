# S3norm
## The second version of the S3norm has been incorporated into the S3V2_IDEAS_ESMP pipeline (https://github.com/guanjue/S3V2_IDEAS_ESMP#Prerequisites-and-S3V2_IDEAS_ESMP-installation). 
### The S3V2_IDEAS_ESMP pipeline is much easier to use and it also incorporate the IDEAS package to call epigenetic state or master peak list across different epigenomic data.

#### From  Reads Count (RC) to S3norm normalized -log10(p-value)

#### Motivation: The quantitative comparison of epigenomic data across multiple cell types has become a promising way to understand the biological function of epigenetic modifications. Due to difference in sequencing depth and signal-to-noise ratio, however, the raw epigenomic data may not reflect the real biological difference between cell types. Existing normalization methods are mainly designed for scaling signals in either the whole-genome or the peak regions, without considering the potentially different scaling factors between peak and background regions. Results: We propose a new data normalization method, S3norm, that normalizes the data by using a monotonic nonlinear data transformation to match signals in both the peak regions and the background regions differently, such that both sequencing depth and signal-to-noise ratio between data sets can be simulatenously normalized. We show that the S3norm normalized epigenomic data can better reflect real biological differences across multiple cell types.



<img src="https://github.com/guanjue/S3norm/blob/master/example_figures/overall_pipeline.png" width="800"/>

##### Figure 1. The overall workflow of the S3norm normalization method. There are three major steps in S3norm. (a) The 1st step is to convert reads count in the 200-bp bins to -log10(p-value) for each epigenomic dataset. Each box represents the different signal tracks of the same data. The first one is the raw reads count of the G1E H3K4me3 dataset. The second one is the reads count of input sample. The third one is the -log10(p-value) of the G1E H3K4me3 dataset. The shoulder of the peaks are reduced after convert reads count to -log10(p-value) with the background adjustment. (b) The 2nd step is selecting the dataset with the highest SNR as the reference dataset for the S3norm normalization. The barplot represents the SNRs of all datasets. The dataset with the highest SNR (dataset with the orange bar) will be selected as the reference dataset. (c) The 3rd step is using a monotonic nonlinear data transformation model to normalize both the SNR and SD between the two datasets. The (1) part is identifing common peak regions and the common background regions between the two datasets. The left scatterplot is showing the signal of each bin in the target dataset and the reference dataset. In the right scatterplot, each data point is colored based on the type of the data point. The orange data points represent the common peak bins. The gray data points represent the common background bins. The blue data points represent the dataset-specific bins. The (2) part is using the monotonic nonlinear data transformation model to rotate the signal of the target dataset, so that (i) the means of the common peak regions of two datasets and (ii) the means of the common background regions of the two datasets can be matched. 

#####################################################################################

## Table of Contents
**[(1) Prerequisites and S3norm installation](#Prerequisites-and-S3norm-installation)**<br>
#####
**[(2) Inputs for S3norm](#Inputs-for-S3norm)**<br>
#####
**[(3) How to run S3norm pipeline](#How-to-run-S3norm-pipeline)**<br>
#####
**[(4) Outputs of S3norm](#Outputs-of-S3norm)**<br>
#####
**[(5) How to run specific steps in S3norm pipeline](#How-to-run-specific-steps-in-S3norm-pipeline)**<br>
#####
**[(6) Contacts and References](#Contacts-and-References)**<br>
#####

#####################################################################################

## Prerequisites and S3norm installation
### S3norm dependencies are as follows:
#### python/2.7 (https://www.python.org/downloads/release/python-2716/)
#### python dependencies: numpy, scipy
#### R (https://www.r-project.org/)
#### gawk

### Installing S3norm pipeline
#### Clone the github repository 
```
cd /where_user_clone_the_S3norm_GitHub/
git clone https://github.com/guanjue/S3norm.git
```
#### Install dependency: 
##### If some of dependencies were not installed, user can use the following command to install them
```
###### For python dependencies, they can be installed by the following scripts
pip install --upgrade pip --user
pip install --upgrade numpy --user
pip install --upgrade scipy --user

###### Installing gwak
### Installing brew
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
### For some MAC, the following script needs to be run before installing gawk by brew
sudo chown -R "$USER":admin $(brew --prefix)/*

### Installing gawk
brew install gawk

```


#####################################################################################

## Inputs for S3norm
#### (1) The input filelist for S3norm
##### The filelist contains the names of the ChIP bedgraph file and the control bedgraph. Each row is one ChIP-seq sample. The 1st column is the ChIP bedgraph and 2nd column is the Contrl bedgraph file. (Separated by tab "\t") 
##### The example of the filelist is in the 'example_file' folder.
##### The first column is the bedgraph file of ChIP signal and second column is the corresponding bedgraph file of the Control signal track.
##### For the ChIP-seq, the Control signal can be computed by the same way in MACS (Zhang, Yong, et al. "Model-based analysis of ChIP-Seq (MACS)." Genome biology 9.9 (2008): R137.)
##### For the ATAC-seq (or any other signal without control), a bedgraph file with control signal all equal to 1 can be used.
```
head file_list.txt 
sig1.sorted.bedgraph	sig1.ctrl.sorted.bedgraph
sig2.sorted.bedgraph	sig2.ctrl.sorted.bedgraph
sig3.sorted.bedgraph	sig3.ctrl.sorted.bedgraph
```

#### (2) The bedgraph files for S3norm
##### The format of bedgraph is as follows:
##### The four columns are chromosome, bin_start, bin_end, and signal in the bin.
##### For the S3norm full pipeline, the average read counts for each bin should be used as signal.
##### For the bedgraph files, they can be generated from bed file AND bigwig files by the bigWigAverageOverBed in UCSC utilities (http://hgdownload.soe.ucsc.edu/admin/exe/)
```
head sig1.UNsorted.bedgraph
chr8	65127400	65127600	77.25
chr21	40481600	40481800	72.84
chr17	19170200	19170400	63.21
chr14	32630400	32630600	51.64
chr6	118552200	118552400	129.82
chr13	93149400	93149600	295.07
chr10	117806400	117806600	142.97
chr19	1370200	1370400	223.43
chr14	28469600	28469800	167.98
chr2	181514400	181514600	220.7
```

#### !!! The first three columns of All of bedgraph files in the filelist should be exactly the same. !!!
#### !!! Only the fourth column is different. !!!
####
#### !!! Each of bedgraph file in the filelist (both ChIP bedgraph and Control bedgraph) need to be sorted before running S3norm. !!!
##### This can be done by the following command:
```
sort -k1,1 -k2,2n sig1.UNsorted.bedgraph > sig1.sorted.bedgraph
sort -k1,1 -k2,2n sig2.UNsorted.bedgraph > sig2.sorted.bedgraph
sort -k1,1 -k2,2n sig3.UNsorted.bedgraph > sig3.sorted.bedgraph
sort -k1,1 -k2,2n sig1.ctrl.UNsorted.bedgraph > sig1.ctrl.sorted.bedgraph
sort -k1,1 -k2,2n sig2.ctrl.UNsorted.bedgraph > sig2.ctrl.sorted.bedgraph
sort -k1,1 -k2,2n sig3.ctrl.UNsorted.bedgraph > sig3.ctrl.sorted.bedgraph

###### The head of the bedgraph files after sorting. 
###### The first three columns of the bedgraph files are exactly the same. 
###### Only the 4th columns are different. 
head sig1.sorted.bedgraph
chr1	7000	7200	0
chr1	18800	19000	0
chr1	62400	62600	5.02
chr1	63800	64000	188.21
chr1	95600	95800	16.41
chr1	136000	136200	0
chr1	156000	156200	0
chr1	158800	159000	0
chr1	206400	206600	51.87
chr1	217000	217200	0

head sig2.sorted.bedgraph
chr1	7000	7200	0
chr1	18800	19000	0
chr1	62400	62600	0
chr1	63800	64000	2.66
chr1	95600	95800	0
chr1	136000	136200	50.26
chr1	156000	156200	0
chr1	158800	159000	0
chr1	206400	206600	0
chr1	217000	217200	0

head sig3.sorted.bedgraph
chr1	7000	7200	0
chr1	18800	19000	0
chr1	62400	62600	0
chr1	63800	64000	0
chr1	95600	95800	0
chr1	136000	136200	0
chr1	156000	156200	0
chr1	158800	159000	0
chr1	206400	206600	0
chr1	217000	217200	0

```


#####################################################################################

## How to run S3norm pipeline
#### Use 's3norm_pipeline.py' to run S3norm pipeline.
##### After perparing the input data, user just need to set the 'script_directory' and 'working_directory' to run S3norm.
##### For 'script_directory', it should be set as the location of folder where the 'S3norm' is saved
##### For 'working_directory', it should be set as the location of folder where the input data (bedgraph files and file_list.txt) are saved.
##### The example script:
```
### Setting script directory
script_directory='/where_user_clone_the_S3norm_GitHub/S3norm/'
### Setting working directory
working_directory='/where_user_clone_the_S3norm_GitHub/S3norm/example_file/'
### Entering working directory
cd $working_directory
### Run S3norm
time python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt
```
#### The same script is also in the 'run_pipeline.sh' in the S3norm folder.
##### To run the 's3norm_pipeline.py' using the 'run_pipeline.sh', user also need to change the 'script_directory' and 'working_directory' in the 'run_pipeline.sh'.
##### Then Run:
```
bash run_pipeline.sh
```
#### The output should looks as follows in the 'working_directory' ('example_file/' folder in the example)
```
ls -ltrh example_file/
total 91464
-rw-r--r--  1 universe  staff   2.8M Jul 29 00:46 sig1.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   3.9M Jul 29 00:47 sig1.ctrl.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   2.7M Jul 29 00:48 sig2.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   3.9M Jul 29 00:48 sig2.ctrl.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   2.7M Jul 29 00:49 sig3.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   2.5M Jul 29 00:49 sig3.ctrl.UNsorted.bedgraph
-rw-r--r--  1 universe  staff   2.8M Jul 29 00:50 sig1.sorted.bedgraph
-rw-r--r--  1 universe  staff   2.7M Jul 29 00:50 sig2.sorted.bedgraph
-rw-r--r--  1 universe  staff   2.7M Jul 29 00:50 sig3.sorted.bedgraph
-rw-r--r--  1 universe  staff   3.9M Jul 29 00:50 sig1.ctrl.sorted.bedgraph
-rw-r--r--  1 universe  staff   3.9M Jul 29 00:50 sig2.ctrl.sorted.bedgraph
-rw-r--r--  1 universe  staff   2.5M Jul 29 00:51 sig3.ctrl.sorted.bedgraph
-rw-r--r--  1 universe  staff   141B Jul 29 00:51 file_list.txt
drwxr-xr-x  4 universe  staff   136B Jul 29 00:52 average_ref_bedgraph
drwxr-xr-x  8 universe  staff   272B Jul 29 00:52 S3norm_rc_bedgraph
drwxr-xr-x  8 universe  staff   272B Jul 29 00:52 S3norm_NBP_bedgraph
drwxr-xr-x  5 universe  staff   170B Jul 29 00:52 NBP_bedgraph
```


#####################################################################################

## Outputs of S3norm
### All outputs will be saved in four subfolders in the working directory.
##### The four subfolders will be named as:
##### 'S3norm_rc_bedgraph/'
##### 'NBP_bedgraph/'
##### 'S3norm_NBP_bedgraph/'
##### 'average_ref_bedgraph/'

#### Within each subfolder, the normalized signals will be saved in '*.begraph' files.
#### The normalization factors will be saved in '*.info.txt' files.
### There are three kinds of outputs for S3norm
##### (1) The S3norm normalized read counts. (Saved in 'S3norm_rc_bedgraph/')
##### (This is the signal for application requires counts data. (e.g. EdgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html) and DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for differential peak calling))
```
ls -ltrh S3norm_rc_bedgraph/
total 22176
-rw-r--r--  1 universe  staff   3.9M Jul 29 00:51 sig1.sorted.bedgraph.s3norm.bedgraph
-rw-r--r--  1 universe  staff    86B Jul 29 00:51 sig1.sorted.bedgraph.info.txt
-rw-r--r--  1 universe  staff   3.5M Jul 29 00:51 sig2.sorted.bedgraph.s3norm.bedgraph
-rw-r--r--  1 universe  staff    86B Jul 29 00:51 sig2.sorted.bedgraph.info.txt
-rw-r--r--  1 universe  staff   3.4M Jul 29 00:51 sig3.sorted.bedgraph.s3norm.bedgraph
-rw-r--r--  1 universe  staff    71B Jul 29 00:51 sig3.sorted.bedgraph.info.txt
```

##### (2) The negative log10 p-value of S3norm normalized read counts based on a negative binomial background model. (Saved in 'NBP_bedgraph/')
##### (This is the signal for peak calling 'bdgpeakcall' and 'bdgbroadcall' in MACS2 (https://github.com/taoliu/MACS))
```
ls -ltrh NBP_bedgraph
total 21480
-rw-r--r--  1 universe  staff   3.8M Jul 29 00:51 sig1.sorted.bedgraph.s3norm.NB.neglog10p.bedgraph
-rw-r--r--  1 universe  staff   3.4M Jul 29 00:52 sig2.sorted.bedgraph.s3norm.NB.neglog10p.bedgraph
-rw-r--r--  1 universe  staff   3.3M Jul 29 00:52 sig3.sorted.bedgraph.s3norm.NB.neglog10p.bedgraph
```

##### (3) The S3norm normalized negative log10 p-value based on a negative binomial background model. (Saved in 'S3norm_NBP_bedgraph/')
##### (This is the signal for genome segmentation (https://github.com/guanjue/IDEAS_2018) 
##### and peak calling by 'bdgpeakcall' and 'bdgbroadcall' in MACS2 (https://github.com/taoliu/MACS))
```
ls -ltrh S3norm_NBP_bedgraph/
total 22480
-rw-r--r--  1 universe  staff   4.0M Jul 29 00:52 sig1.sorted.bedgraph.NBP.s3norm.bedgraph
-rw-r--r--  1 universe  staff    71B Jul 29 00:52 sig1.sorted.bedgraph.NBP.info.txt
-rw-r--r--  1 universe  staff   3.6M Jul 29 00:52 sig2.sorted.bedgraph.NBP.s3norm.bedgraph
-rw-r--r--  1 universe  staff    86B Jul 29 00:52 sig2.sorted.bedgraph.NBP.info.txt
-rw-r--r--  1 universe  staff   3.4M Jul 29 00:52 sig3.sorted.bedgraph.NBP.s3norm.bedgraph
-rw-r--r--  1 universe  staff    86B Jul 29 00:52 sig3.sorted.bedgraph.NBP.info.txt
```
##### !!! Noted that for genome segmentation, we suggest user to user use '-r median' for computing reference track.
##### !!! This parameter can be set in last line of the "run_pipeline.sh" script
```
time python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt -r median
```
##### (4) The 4th folder is used to save the reference signals for S3norm. (Saved in 'S3norm_NBP_bedgraph/')
```
ls -ltrh average_ref_bedgraph/
total 13296
-rw-r--r--  1 universe  staff   2.7M Jul 29 00:51 average_ref.bedgraph
-rw-r--r--  1 universe  staff   3.8M Jul 29 00:52 average_ref.bedgraph.NBP.bedgraph
``` 


#####################################################################################

### Parameters for S3norm
#### Required parameters:
##### For S3norm, there are just two required parameters.
##### (1) -s script_folder: the script directory of S3norm (e.g. /Users/universe/Documents/2018_BG/S3norm/src/)
##### (2) -s input_file_list: the input filelist for S3norm (e.g. The 'file_list.txt' in the 'S3norm/example_file/' folder)

#### For other parameters, user can change them by using the following command:
```
script_directory='/Users/universe/Documents/2018_BG/S3norm/'
python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt -r max1 -m non0mean -i 2.0 -f 0.05 -l 0.001 -a 100000 -b 0 -p z -k 0 -g 0
```

#### The other parameters that can be changed:
```
(1) -r : The method for generating the reference signal track. Options: max1 (default: select the dataset with the max FRiP score as reference), max1 (select the dataset with the highest FRiP score as reference), median (generate signal track by using the median signal of each bin), mean (generate signal track by using the mean signal of each bin)
(2) -m : The method for matching peaks and background. Options: non0mean (default), non0median, mean, median) 
(2) -m : filelist_row_number (If user want to select one sample as the reference, use "-m filelist_row_number", where filelist_row_number is the row number (start from 1) of the reference sample in the file_list.txt)
(3) -i : The initial value for the power parameter in the non-linear transformation. Default: 2.0
(4) -f : The FDR threshold for identifying common peaks. Default: 0.05 . The range is 0.0 < x < 1.0
(5) -l : The minimum proportion of bins are used as peak for S3norm. Default: 0.001 . The range is 0.0 < x < 1.0
(6) -a : The upperlimit for signal. This is to reduce the bias cause by extrame signals in the data. Default: 100000
(7) -b : The lowerlimit for signal. S3norm requires all signal to be x >= 0 . Default: 0
(8) -p : The method used to identify common peaks. Options: z (Default) and neglog10p (negative log10 p-value from background model)
(9) -k : The user given common peak regions. Options: 0 (Default, the common peak will be identified by S3norm) or filename (a file points out which bins are the common peaks. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common peak, the column should be 1 for that row. Otherwise, it should be 0 )
(10) -g : The user given common background regions. Options: 0 (Default, the common background will be identified by S3norm) or filename (a file points out which bins are the common background. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common background, the column should be 1 for that row. Otherwise, it should be 0 )
(11) -c : Whether to use cross feature mode. (T: for use cross mark mode; F: for NOT use cross mark mode)
```


#####################################################################################
## How to run specific steps in S3norm pipeline
### The S3norm pipeline has two steps can be run separately.
#### (1) Get S3norm normalized read counts
##### There are three required parameters. 
##### -r : The filename of the reference signal bedgraph file. (The signal track S3norm normalize other files TO)
##### -t : The filename of the traget signal bedgraph file. (The signal to be normalized by S3norm)
##### -o : The output filename of the S3norm normalized signal track (bedgraph format)
```
time python $script_directory'/src/s3norm.py' -r $working_directory'average_ref_bedgraph/average_ref.bedgraph' -t sig1.sort.bedgraph -o sig1.output
```
##### The other parameters can be changed in 's3norm.py' are as follows
```
time python $script_directory'/src/s3norm.py' -r $working_directory'average_ref_bedgraph/average_ref.bedgraph' -t sig1.sorted.bedgraph -o sig1.runseparately.output -m non0mean -i 2.0 -f 0.05 -l 0.001 -a 100000 -b 0 -p z -k 0 -g 0 -c F
```
```
(1) -m : The method for matching peaks and background. Options: non0mean (default), non0median, mean, median)
(2) -i : The initial value for the power parameter in the non-linear transformation. Default: 2.0
(3) -f : The FDR threshold for identifying common peaks. Default: 0.05 . The range is 0.0 < x < 1.0
(4) -l : The minimum proportion of bins are used as peak for S3norm. Default: 0.001 . The range is 0.0 < x < 1.0
(5) -a : The upperlimit for signal. This is to reduce the bias cause by extrame signals in the data. Default: 100000
(6) -b : The lowerlimit for signal. S3norm requires all signal to be x >= 0 . Default: 0
(7) -p : The method used to identify common peaks. Options: z (Default) and neglog10p (negative log10 p-value from background model)
(8) -k : The user given common peak regions. Options: 0 (Default, the common peak will be identified by S3norm) or filename (a file points out which bins are the common peaks. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common peak, the column should be 1 for that row. Otherwise, it should be 0 )
(9) -g  : The user given common background regions. Options: 0 (Default, the common background will be identified by S3norm) or filename (a file points out which bins are the common background. The rows in this file match the rows in bedgraph files. It should contain only 1 column. If the row of a bin is a common background, the column should be 1 for that row. Otherwise, it should be 0 )
(10) -c : Whether to use cross feature mode. (T: for use cross mark mode; F: for not use cross mark mode)
```


#### (2) Get signal track of negative log10 p-value based on a NB background model (NBP)
##### There are three required parameters separated by white space. 
##### 1st: The filename of ChIP signal track after S3nom (bedgraph format).
##### 2nd: The filename of Control signal track after S3nom (bedgraph format).
##### 3rd: The output filename of the NBP signal track (bedgraph format).
```
Rscript $script_directory'/src/negative_binomial_neglog10p.R' $working_directory'S3norm_rc_bedgraph/sig1.sorted.bedgraph.s3norm.bedgraph' sig1.ctrl.sorted.bedgraph sig1.sorted.bedgraph.s3norm.NB.neglog10p.bedgraph
```



#####################################################################################

## Contacts and References
#### Contacts: 
##### gzx103@psu.edu

#### S3norm
Xiang, Guanjue, et al. "S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data." bioRxiv (2018): 506634.

#### Vision paper

