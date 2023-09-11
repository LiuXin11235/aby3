
log_folder_list=(./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_average ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_cipher_index ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_new_search ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_max ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric ./Record/Record_metric)
task_list=("average" "average" "average" "average" "cipher_index" "cipher_index" "cipher_index" "cipher_index" "new_search" "new_search" "new_search" "new_search" "max" "max" "max" "max" "metric" "metric" "metric" "metric")
M_list=(1 1 1 1 1 1 1 1 1 1 1 1 1024 16384 32768 65536 1 1 1 1)
N_list=(1048576 16777216 134217728 1073741824 1048576 16777216 134217728 1073741824 1048576 16777216 134217728 1073741824 1024 16384 32768 65536 1048576 16777216 134217728 1073741824)
c_list=(8 39 256 256 24 112 256 256 27 115 256 256 24 256 256 256 51 136 256 256 )
optB_list=(131072 360448 360448 360448 43690 149796 524288 737280 38836 145888 524288 892928 43690 647168 647168 647168 20560 123361 524288 2247671)


test_times=2;outLimit=10;retry_threshold=5;K=1;

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
    done;
  done;
done;