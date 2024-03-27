# compile the main.
task=$1
args_list=$2

echo ${args_list}

./Eval/dis_mpi_exec.sh ${task_num} "${test_args}"