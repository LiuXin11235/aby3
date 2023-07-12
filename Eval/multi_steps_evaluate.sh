# task_list=("sort" "max" "min" "medium" "rank")
# log_folder_list=(./Record/Record_sort ./Record/Record_max ./Record/Record_min ./Record/Record_medium ./Record/Record_rank)
task_list=("rank" "sort")
log_folder_list=(./Record/Record_rank ./Record/Record_sort)
# task_list=("sort")
# log_folder_list=(./Record/Record_sort)

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

# N_list=(1024 16384 32768 65536)
N_list=(32768)
M=1; K=1
# task_num_list=(256 128 64 32 16 4 1)
task_num_list=(64 128 256)
optB_list=(16388 65536 262144 1048576 4194304)
# optB_list=(65536 262144 1048576 4194304)
exceed_time=(15 10 8 6 3)
repeat=1; test_times=10; retry_threshold=10;


for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for (( t=0; t<${#task_num_list[@]}; t++ )); do
    taskN=${task_num_list[t]}; outLimit=${exceed_time[t]};
    for optB in ${optB_list[@]}; do
      for (( k=0; k<$test_times; k++ )); do
        for N in ${N_list[@]}; do
          # run the tasks with retrying.
          j=0;
          while [ $j -lt $retry_threshold ]; do
            echo "retrying? "$j" parameters: taskN="$taskN" n="$N" m="$M" k="$K" repeat="$repeat" task="$task" optB="$optB >> ${log_folder}/record.log;
            timeout ${outLimit}m ./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M k=$K repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
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
    done;
  done;
  # tar the results.
  tar cvf $log_folder.tar $log_folder;
done;
