task_list=("average" "average" "average" "average" "cipher_index" "cipher_index" "cipher_index" "cipher_index" "new_search" "new_search" "new_search" "new_search" "max" "max" "max" "max" "metric" "metric" "metric" "metric")
log_folder_list=(./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric)
M_list=(1 1 1 1 1 1 1 1 1 1 1 1 1024 16384 32768 65536 1 1 1 1)
N_list=(1048576 16777216 134217728 1073741824 1048576 16777216 134217728 1073741824 1048576 16777216 134217728 1073741824 1024 16384 32768 65536 1048576 16777216 134217728 1073741824)
c_list=(2 28 256 256 15 100 256 256 19 111 256 256 14 256 256 256 57 127 256 256)
optB_list=(360448 360448 360448 360448 69905 167772 221184 221184 55188 122880 122880 122880 74898 647168 647168 647168 18396 132104 524288 3858432)

test_times=2;exceed_time=10;retry_threshold=5;K=1;

for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done


for (( i=0; i<${#task_list[@]}; i++ )); do
  task=${task_list[i]}; log_folder=${log_folder_list[i]}; 
  N=${N_list[i]}; M=${M_list[i]};
  taskN=${c_list[i]}; optB=${optB_list[i]};

  for (( z=0; z<${test_times}; z++ )); do
    j=0;
    while [ $j -lt $retry_threshold ]; do
      timeout ${outLimit}m ./Eval/mpi_subtask.sh taskN=${taskN} n=$N m=$M k=$K task=$task optB=$optB log_folder=$log_folder/;
      if [ $? -eq 0 ]; then
        break; 
      fi
      ./Eval/kill_all.sh frontend;
      j=$(expr $j + 1);
      if [ $j -eq $retry_threshold ]; then
        echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
      fi
  done
done
