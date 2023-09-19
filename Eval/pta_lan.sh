record_list=(./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_max)
task_list=(new_search new_search new_search new_search new_search new_search metric metric metric metric metric metric average average average average average average cipher_index cipher_index cipher_index cipher_index cipher_index cipher_index max max max max max max)
n_list=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1024 2048 4096 8192 16384 32768)
m_list=(1048576 4194304 16777216 67108864 268435456 1073741824 1048576 4194304 16777216 67108864 268435456 1073741824 1048576 4194304 16777216 67108864 268435456 1073741824 1048576 4194304 16777216 67108864 268435456 1073741824 1024 2048 4096 8192 16384 32768)
c_list=(27 54 115 235 256 256 51 51 136 256 256 256 8 19 39 79 256 256 24 53 112 229 256 256 24 60 122 248 256 256)
B_list=(38837 77673 145889 285570 892928 892928 20561 82242 123362 262144 1048576 2247671 131072 220753 360448 360448 360448 360448 43691 79138 149797 293052 737280 737280 43691 69906 137519 270601 647168 647168)

test_times=10;outLimit=10;retry_threshold=5;K=1;

./Eval/basic/seq_test.sh \
  "$(IFS=":"; echo "${task_list[*]}")" \
  "$(IFS=":"; echo "${record_list[*]}")" \
  "$(IFS=":"; echo "${n_list[*]}")" \
  "$(IFS=":"; echo "${m_list[*]}")" \
  "$(IFS=":"; echo "${c_list[*]}")" \
  "$(IFS=":"; echo "${B_list[*]}")" \
  ${test_times} ${outLimit} ${retry_threshold} ${K}