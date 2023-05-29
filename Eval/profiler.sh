task="vector"
log_folder=./Record/Record_vector

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# backup the last log folder.
if [ ! -d ${log_folder} ]; then
  mkdir ${log_folder}
else
  cp -r ${log_folder} ${log_folder}${day}${timeStamp};
  rm -r ${log_folder};
  mkdir ${log_folder};
fi

# compile
python build.py

# synchroonize with others
scp ./bin/frontend aby31:~/aby3/bin/ &
scp ./bin/frontend aby32:~/aby3/bin/ &
wait;

tasks_list=(8 16 32 64 128)
# n_list=(10 100 1000 10000 100000 1000000 10000000)
# repeats_list=(1000 500 500 50 10 5 1)

for taskN in ${tasks_list[@]}; do

  day=$(date +%m-%d);
  timeStamp=$(date +"%H%M%s");
  echo ${day}-${timeStamp}-${taskN}-${n} >> ${log_folder}/log.time

  mpirun -n ${taskN} ./bin/frontend -prog -1 -role 0 -N 0 -TASK_NUM ${taskN} -FUNC ${task} -M 0 -OPT_BLOCK -1 -repeats 0 &
  ssh aby31 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 1 -N 0 -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats 0 " &
  ssh aby32 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 2 -N 0 -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats 0 " &
  wait;

  ./Eval/kill_all frontend;
  sleep 5;

done;

# for taskN in ${tasks_list[@]}; do
#   for (( i=0; i<${#n_list[@]}; i++ )); do
#     n=${n_list[i]};
#     repeat=${repeats_list[i]};
    
#     day=$(date +%m-%d);
#     timeStamp=$(date +"%H%M%s");
#     echo ${day}-${timeStamp}-${taskN}-${n} >> ${log_folder}/log.time

#     mpirun -n ${taskN} ./bin/frontend -prog -1 -role 0 -N ${n} -TASK_NUM ${taskN} -FUNC ${task} -M 0 -OPT_BLOCK -1 -repeats ${repeat} &
#     ssh aby31 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 1 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" >> ./log" &
#     ssh aby32 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 2 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" >> ./log" &
#     wait;

#     sleep 60;

#   done;
# done;

