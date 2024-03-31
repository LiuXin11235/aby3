task_num=$1
args_list=$2
echo ${args_list}

USER_FOLDER=/home/tsingj_ubuntu/fanxy/PtA/aby3

ulimit -n 65536;
mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 0 ${args_list} >> ./log 2>&1 &
# mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 1 ${args_list} &
# mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 2 ${args_list} &
ssh sosp1 "ulimit -n 65536; mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 1 ${args_list} &" &
ssh sosp2 "ulimit -n 65536; mpirun -np $task_num $USER_FOLDER/out/build/linux/frontend/frontend -role 2 ${args_list} &" &

wait;