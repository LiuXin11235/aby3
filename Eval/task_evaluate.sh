task_list=("index" "rank" "search")
log_folder_list=(./Record/Record_index ./Record/Record_rank ./Record/Record_search)

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
N=1000
M=1
TASKS=64
OPT_B=500

for task in ${task_list[@]}; do

  echo ${task}
  mpirun -n ${TASKS} ./bin/frontend -prog -1 -role 0 -testFlag ${testFlag} -N ${N} -M ${M} -TASK_NUM ${TASKS} -OPT_BLOCK ${OPT_B} -FUNC ${task} >> ./log &

  ssh aby31 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog -1 -role 1 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" -FUNC ${task} >> ./log" &

  ssh aby32 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog -1 -role 2 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" -FUNC ${task} >> ./log" &

  wait;

done;
