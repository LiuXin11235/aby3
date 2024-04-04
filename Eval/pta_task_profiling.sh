task_list=("cipher_index" "max" "sort" "metric" "sum")
# task_list=("cipher_index")

for task in ${task_list[@]}; do
    python ./PtA_deploy/system_profile.py --task_profile ${task}
done