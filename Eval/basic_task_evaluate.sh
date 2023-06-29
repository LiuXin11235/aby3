# 1d functions.
# task_list=("average" "index" "search" "new_search" "select")
task_list=("average")
log_folder_list=(./Record/Record_average)
# task_list=("select")
# log_folder_list=(./Record/Record_select)
# log_folder_list=(./Record/Record_average ./Record/Record_index ./Record/Record_search ./Record/Record_new_search ./Record/Record_select)


day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# cleanup the old folder and create a clean one.
for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done

# compile
python build.py

# synchroonize with others
scp ./bin/frontend aby31:~/aby3/bin &
scp ./bin/frontend aby32:~/aby3/bin &
wait;

# test settings.
# N_list=(1048576 16777216 134217728 1073741824)
N_list=(1073741824)
# N_list=(500000000)
repeat=3; test_times=3; retry_threshold=5
task_num_list=(256 128)
# task_num_list=(256 128 64 32 16 4 1)
optB_list=(1024 16384 131072 1048576 16777216 134217728)
# optB_list=(1048576 1048576 1048576 1048576 1048576 1048576 1048576)
exceed_time=(50 100)

for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for (( t=0; t<${#task_num_list[@]}; t++ )); do
    taskN=${task_num_list[t]}; outLimit=${exceed_time[t]};
    M=1; K=4
    for optB in ${optB_list[@]}; do
      for (( k=0; k<$test_times; k++ )); do
        for N in ${N_list[@]}; do
          # run the tasks with retrying.
          j=0;
          while [ $j -lt $retry_threshold ]; do
            timeout ${outLimit}m ./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M k=$K repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
            if [ $? -eq 0 ]; then
              break; 
            fi
            j=$(expr $j + 1);
            if [ $j -eq $retry_threshold ]; then
              echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
            fi
          done;
        done;
      done;
    done;
  done;
  tar cvf $log_folder.tar $log_folder;
done;