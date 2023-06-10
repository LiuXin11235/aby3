task_list=("search" "select" "index" "average")
log_folder_list=(./Record/Record_search ./Record/Record_select ./Record/Record_index ./Record/Record_average)

# task_list=("average")
# log_folder_list=(./Record/Record_average)

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
scp ./bin/frontend aby31:~/aby3/bin/ &
scp ./bin/frontend aby32:~/aby3/bin/ &
wait;

# test settings.
N_list=(1000000)
M=1
optB_list=(250000)
repeat=1; test_times=10; retry_threshold=5
task_num_list=(16)


for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for taskN in ${task_num_list[@]}; do
    for (( k=0; k<$test_times; k++ )); do
      for N in ${N_list[@]}; do
        for optB in ${optB_list[@]}; do
          # run the tasks with retrying.
          j=0;
          while [ $j -lt $retry_threshold ]; do
            timeout 15m ./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
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

  # tar the results.
  # tar cvf $log_folder.tar $log_folder;
done;