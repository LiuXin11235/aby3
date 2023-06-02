task_list=("index" "search")
log_folder_list=(./Record/Record_index ./Record/Record_search)

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
N=1000000000
M=1
optB=80000
repeat=1

task_num_list=(256 180)

for (( i=0; i<${#task_list[@]}; i++ )); do
  for taskN in ${task_num_list[@]}; do

    task=${task_list[i]};
    log_folder=${log_folder_list[i]};

    j=0;
    while [ $j -lt 5 ]; do
      timeout 50m ./Eval/general_mpi.sh ${taskN} ${N} ${M} ${repeat} ${task} ${optB} ${log_folder}/
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
    # mpirun -n ${TASKS} ./bin/frontend -prog -1 -role 0 -testFlag ${testFlag} -N ${N} -M ${M} -TASK_NUM ${TASKS} -OPT_BLOCK ${OPT_B} -FUNC ${task} >> ./log &

    # ssh aby31 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog -1 -role 1 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" -FUNC ${task} >> ./log" &

    # ssh aby32 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog -1 -role 2 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" -FUNC ${task} >> ./log" &
    wait;
  done;
done;
