log_folder_list=(./Record/Record_average ./Record/Record_cipher_index ./Record/Record_new_search ./Record/Record_max ./Record/Record_metric)
task_list=("average" "cipher_index" "new_search" "max" "metric")
M_list=(1 1 1 1024 1)
N_list=(65536 65536 65536 1024 65536)
c_list=(64 64 64 64 64)
optB_list=(256 256 256 256 256)

test_times=1;outLimit=10;retry_threshold=5;K=1;

./Eval/basic/seq_test.sh \
  "$(IFS=":"; echo "${task_list[*]}")" \
  "$(IFS=":"; echo "${log_folder_list[*]}")" \
  "$(IFS=":"; echo "${M_list[*]}")" \
  "$(IFS=":"; echo "${N_list[*]}")" \
  "$(IFS=":"; echo "${c_list[*]}")" \
  "$(IFS=":"; echo "${optB_list[*]}")" \
  ${test_times} ${outLimit} ${retry_threshold} ${K}