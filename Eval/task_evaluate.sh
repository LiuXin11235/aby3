# task_list=("search" "select" "index")
# log_folder_list=(./Record/Record_search ./Record/Record_select ./Record/Record_index)
task_list=("search")
log_folder_list=(./Record/Record_search)

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# backup the old folder and creat a clean folder.
for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    cp -r ${log_folder} ${log_folder}${day}${timeStamp};
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

#test
N_list=(1000000000)
M=1
optB_list=(50000 100000 500000 1000000)
repeat=1
testing_times=(1 1 1 1 1)

task_num_list=(256)

for test in ${testing_times[@]}; do

  for N in ${N_list[@]}; do
    for optB in ${optB_list[@]}; do
      for (( i=0; i<${#task_list[@]}; i++ )); do
        for taskN in ${task_num_list[@]}; do

          task=${task_list[i]};
          log_folder=${log_folder_list[i]};

          j=0;
          while [ $j -lt 5 ]; do
            timeout 15m ./Eval/general_mpi.sh ${taskN} ${N} ${M} ${repeat} ${task} ${optB} ${log_folder}/
            if [ $? -eq 0 ];then
              break;
            fi
            echo ">>>>>>>> in retry!"${taskN}-${optB};
            ./Eval/kill_all.sh frontend;
            sleep 1m;
            echo ">>>>>>>> after sleep!"${taskN}-${optB};
            ipcs -m shmid;
            echo "retry: "${j}" for "${optB}-${taskN} >> ${log_folder}/error.log;

            rm ${log_folder}/log-config-N=${n}-TASKS=${taskN}-Vec=${optB}-0;

            j=$(expr $j + 1);
            if [ $j -eq 5 ];then
              echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
            fi
          done;
          wait;
        done;
      done;
    done;
  done;
done;


# task_list=("rank")
# log_folder_list=(./Record/Record_rank)

# day=$(date +%m-%d);
# timeStamp=$(date +"%H%M%s");

# # backup the old folder and creat a clean folder.
# for log_folder in ${log_folder_list[@]}; do
#   if [ ! -d ${log_folder} ]; then
#     mkdir ${log_folder}
#   else
#     cp -r ${log_folder} ${log_folder}${day}${timeStamp};
#     rm -r ${log_folder};
#     mkdir ${log_folder};
#   fi
# done

# # compile
# python build.py

# # synchroonize with others
# scp ./bin/frontend aby31:~/aby3/bin/ &
# scp ./bin/frontend aby32:~/aby3/bin/ &
# wait;

# #test
# N_list=(100000)
# M=1
# optB_list=(5000 10000 50000 100000 500000 1000000 5000000 10000000)
# repeat=1
# testing_times=(1 1 1 1 1)

# task_num_list=(256)

# for test in ${testing_times[@]}; do

#   for N in ${N_list[@]}; do
#     for optB in ${optB_list[@]}; do
#       for (( i=0; i<${#task_list[@]}; i++ )); do
#         for taskN in ${task_num_list[@]}; do

#           task=${task_list[i]};
#           log_folder=${log_folder_list[i]};

#           j=0;
#           while [ $j -lt 5 ]; do
#             timeout 100m ./Eval/general_mpi.sh ${taskN} ${N} ${N} ${repeat} ${task} ${optB} ${log_folder}/
#             if [ $? -eq 0 ];then
#               break;
#             fi
#             echo ">>>>>>>> in retry!"${taskN}-${optB};
#             ./Eval/kill_all.sh frontend;
#             sleep 1m;
#             echo ">>>>>>>> after sleep!"${taskN}-${optB};
#             ipcs -m shmid;
#             echo "retry: "${j}" for "${optB}-${taskN} >> ${log_folder}/error.log;

#             rm ${log_folder}/log-config-N=${n}-TASKS=${taskN}-Vec=${optB}-0;

#             j=$(expr $j + 1);
#             if [ $j -eq 5 ];then
#               echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
#             fi
#           done;
#           wait;
#         done;
#       done;
#     done;
#   done;
# done;
