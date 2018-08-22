###### set parameters
script_dir=/storage/home/gzx103/group/software/vision_mouse/src/
working_dir=/storage/home/gzx103/scratch/vision/test_pipeline/
input_dir=/storage/home/gzx103/scratch/vision/test_pipeline/input_5end_rc/
input_file_list=info_table_all.rc2nbp.txt
overall_upper=16
overall_lower=2
select_method=frip
user_given_global_ref=ERY_fl.h3k4me3rep.fisher_p.txt

time bash $script_dir'overall_pipeline.sh' $script_dir $working_dir $input_dir $input_file_list $overall_upper $overall_lower $select_method $user_given_global_ref

