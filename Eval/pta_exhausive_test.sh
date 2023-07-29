# task_list=("metric" "new_search")
# log_folder_list=(./Record/Record_metric_all ./Record/Record_new_search_all)

# # cleanup the old folder and create a clean one.
# for log_folder in ${log_folder_list[@]}; do
#   if [ ! -d ${log_folder} ]; then
#     mkdir ${log_folder}
#   else
#     rm -r ${log_folder};
#     mkdir ${log_folder};
#   fi
# done

# # compile
# python build.py

# # synchroonize with others
# scp ./bin/frontend aby31:~/aby3/bin &
# scp ./bin/frontend aby32:~/aby3/bin &
# wait;

test_times=2;

# # inputs.
# N_list=(1048576 16777216 134217728 1073741824)
# # N_list=(1073741824)
# M_list=(1 1 1 1)
# K=1

# test configs.
c_list=(256 64 16 4 1)
# c_list=(256)
# optB_list=(65536)
optB_list=(65536 262144 1048576 4194304)
exceed_time_list=(7 10 13 20 80)
retry_threshold=5;


# for (( i=0; i<${#task_list[@]}; i++ )); do

#   #   # set the task the correspondng log folder
#   task=${task_list[i]}; log_folder=${log_folder_list[i]};

#   for (( t=0; t<${#N_list[@]}; t++ )); do
#     N=${N_list[t]}; M=${M_list[t]};

#     for (( c=0; c<${#c_list[@]}; c++ )); do
#       taskN=${c_list[c]}; outLimit=${exceed_time_list[c]};

#       for (( k=0; k<${#optB_list[@]}; k++ )); do
#         optB=${optB_list[k]};

#         for (( z=0; z<${test_times}; z++  )); do
#           j=0;
#           while [ $j -lt $retry_threshold ]; do
#             timeout ${outLimit}m ./Eval/mpi_subtask.sh taskN=${taskN} n=$N m=$M k=$K task=$task optB=$optB log_folder=$log_folder/;
#             if [ $? -eq 0 ]; then
#               break; 
#             fi
#             ./Eval/kill_all.sh frontend;
#             j=$(expr $j + 1);
#             if [ $j -eq $retry_threshold ]; then
#               echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
#             fi
#           done;
#         done;
#       done;
#     done;
#   done;
# done;


# squared task
task_list=("max")
log_folder_list=(./Record/Record_max_all)

for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done

#test_settings.
N_list=(1024 16384 32768 65536)
M_list=(1 1 1 1)
K=1
exceed_time_list=(7 10 13 20 80)
retry_threshold=5;


for (( i=0; i<${#task_list[@]}; i++ )); do

  #   # set the task the correspondng log folder
  task=${task_list[i]}; log_folder=${log_folder_list[i]};

  for (( t=0; t<${#N_list[@]}; t++ )); do
    N=${N_list[t]}; M=${M_list[t]};

    for (( c=0; c<${#c_list[@]}; c++ )); do
      taskN=${c_list[c]}; outLimit=${exceed_time_list[c]};

      for (( k=0; k<${#optB_list[@]}; k++ )); do
        optB=${optB_list[k]};

        for (( z=0; z<${test_times}; z++  )); do
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
    done;
  done;
done;