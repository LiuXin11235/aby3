task_num=$1
args_list=$2
echo ${args_list}

mpirun -np $task_num ./out/build/linux/frontend/frontend -role 0 ${args_list} &
mpirun -np $task_num ./out/build/linux/frontend/frontend -role 1 ${args_list} &
mpirun -np $task_num ./out/build/linux/frontend/frontend -role 2 ${args_list} &
wait;