# task_list=("average" "index" "new_search" "metric")
# log_folder_list=(./Record/Record_average_tgl ./Record/Record_index_tgl ./Record/Record_new_search_tgl ./Record/Record_metric_tgl)

# task_list=("metric")
# log_folder_list=(./Record/Record_metric_tgl)

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

# #test_settings.
# N_list=(1048576 16777216 134217728 1073741824)
# M_list=(1 1 1 1)
# K=1

# for (( i=0; i<${#task_list[@]}; i++ )); do

#   #   # set the task the correspondng log folder
#   task=${task_list[i]};
#   log_folder=${log_folder_list[i]};

#   for (( t=0; t<${#N_list[@]}; t++ )); do
#     N=${N_list[t]}; M=${M_list[t]};
#     taskN=1; optB=$((N * M));
#     ./Eval/mpi_subtask.sh taskN=${taskN} n=$N m=$M k=$K task=$task optB=$optB log_folder=$log_folder/;
#   done
# done


# squared task
task_list=("max")
log_folder_list=(./Record/Record_max_tgl)

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

for (( i=0; i<${#task_list[@]}; i++ )); do

  #   # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for (( t=0; t<${#N_list[@]}; t++ )); do
    N=${N_list[t]}; M=${M_list[t]};
    taskN=1; optB=$((N * N));
    ./Eval/mpi_subtask.sh taskN=${taskN} n=$N m=$M k=$K task=$task optB=$optB log_folder=$log_folder/;
  done
done