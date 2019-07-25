### Setting script directory !!! user needs to change the directory where S3norm installed !!!
script_directory='/Users/universe/Documents/2018_BG/S3norm/'

### Setting working directory !!! user needs to change the directory where input data is saved !!!
working_script_directory='/Users/universe/Documents/2018_BG/S3norm/example_file/'

### Entering working directory
cd $working_script_directory

### Run S3norm
time python $script_directory'/src/S3norm_pipeline.py' -s $script_directory'/src/' -t file_list.txt