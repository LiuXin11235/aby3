# task_list=("bio_metric" "mean_distance")
# log_folder_list=(./Record/Record_bio_metric ./Record/Record_mean_distance)
task_list=("metric")
log_folder_list=(./Record/Record_metric)

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

# N_list=(1048576 16777216 134217728 1073741824)
N_list=(1073741824 134217728 16777216 1048576)
M=1; K=2;
# M_list=(1)
# K_list=(2)
repeat=1; test_times=3; retry_threshold=5
# task_num_list=(256 128 64 32 16 4 1)
task_num_list=(1 4 16 32 64 128 256)
optB_list=(524288 524288 524288 524288 524288 524288 524288)
exceed_time=(100 75 45 30 15 15 15)
# exceed_time=(15 15 15 30 45 75 100)

for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for (( t=0; t<${#task_num_list[@]}; t++ )); do
    taskN=${task_num_list[t]}; optB=${optB_list[t]}; outLimit=${exceed_time[t]};
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
  # tar the results.
  tar cvf $log_folder.tar $log_folder;
done;

