taskN=$1
n=$2
repeat=$3
task="vector"

mpirun -n ${taskN} ./bin/frontend -prog -1 -role 0 -N ${n} -TASK_NUM ${taskN} -FUNC ${task} -M 0 -OPT_BLOCK -1 -repeats ${repeat} &
ssh aby31 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 1 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" " &
ssh aby32 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 2 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" " &
wait;