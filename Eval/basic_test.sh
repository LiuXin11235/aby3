log_folder_list=(./Record/Record_new_search)
task_list=("new_search")
M_list=(1)
N_list=(1073741824)
c_list=(256)
optB_list=(72025)

# log_folder_list=(./Record/Record_metric)
# task_list=("metric")
# M_list=(1)
# N_list=(1073741824)
# c_list=(256)
# optB_list=(4194304)

rm -r ./DEBUG/*

test_times=1;outLimit=10;retry_threshold=5;K=1;

./Eval/basic/seq_test.sh \
  "$(IFS=":"; echo "${task_list[*]}")" \
  "$(IFS=":"; echo "${log_folder_list[*]}")" \
  "$(IFS=":"; echo "${M_list[*]}")" \
  "$(IFS=":"; echo "${N_list[*]}")" \
  "$(IFS=":"; echo "${c_list[*]}")" \
  "$(IFS=":"; echo "${optB_list[*]}")" \
  ${test_times} ${outLimit} ${retry_threshold} ${K}