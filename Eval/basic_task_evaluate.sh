task_list=("index" "bio_metric" "average"  "search" "select" "mean_distance")
log_folder_list=(./Record/Record_index ./Record/Record_bio_metric ./Record/Record_average ./Record/Record_search ./Record/Record_select ./Record/Record_mean_distance)


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
N_list=(268435456 67108864 16777216 4194304 1048576)
repeat=3; test_times=3; retry_threshold=5
task_num_list=(256 128 64 32 16 4)
optB_list=(1048576 1048576 1048576 1048576 1048576)

for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for (( t=0; t<${#task_num_list[@]}; t++ )); do
    optB=${optB_list[t]}
    M=16
    K=4
    taskN=${task_num_list[t]}
    for (( k=0; k<$test_times; k++ )); do
      for N in ${N_list[@]}; do
        # run the tasks with retrying.
        j=0;
        while [ $j -lt $retry_threshold ]; do
          timeout 20m ./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M k=$K repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
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
  tar cvf $log_folder.ptr.tar $log_folder;
done;


# # test settings.
# N_list=(1000000000 100000000 10000000 1000000)
# M_list=(1)
# optB_list=(262144)
# repeat=1; test_times=5; retry_threshold=5
# # task_num_list=(256 128 64 32 16)
# task_num_list=(1)
# # task_num_list=(64)

# # cleanup the old folder and create a clean one.
# for log_folder in ${log_folder_list[@]}; do
#   if [ ! -d ${log_folder} ]; then
#     mkdir ${log_folder}
#   else
#     rm -r ${log_folder};
#     mkdir ${log_folder};
#   fi
# done

# for (( i=0; i<${#task_list[@]}; i++ )); do

#   # set the task the correspondng log folder
#   task=${task_list[i]};
#   log_folder=${log_folder_list[i]};

#   for taskN in ${task_num_list[@]}; do
#     for (( k=0; k<$test_times; k++ )); do
#       for N in ${N_list[@]}; do
#         for M in ${M_list[@]}; do
#         # for M in ${M_list[@]}; do
#           for optB in ${optB_list[@]}; do
#             # run the tasks with retrying.
#             j=0;
#             while [ $j -lt $retry_threshold ]; do
#               timeout 150m ./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
#               if [ $? -eq 0 ]; then
#                 break; 
#               fi
#               j=$(expr $j + 1);
#               if [ $j -eq $retry_threshold ]; then
#                 echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
#               fi
#             done;
#           done;
#         done;
#       done;
#     done;
#   done;

#   # tar the results.
#   tar cvf $log_folder.single.tar $log_folder;
# done;