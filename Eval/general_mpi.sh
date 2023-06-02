taskN=$1
n=$2
m=$3
repeat=$4
task=$5
optB=$6
logFolder=$7

mpirun -bind-to hwthread -np ${taskN} ./bin/frontend -prog -1 -role 0 -N ${n} -TASK_NUM ${taskN} -FUNC ${task} -M ${m} -OPT_BLOCK ${optB} -repeats ${repeat} -logFolder ${logFolder}&
ssh aby31 "cd ./aby3/; mpirun -bind-to hwthread -np "${taskN}" ./bin/frontend -prog -1 -role 1 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M "${m}" -OPT_BLOCK "${optB}" -repeats "${repeat}" -logFolder "${logFolder}" " &
ssh aby32 "cd ./aby3/; mpirun -bind-to hwthread -np "${taskN}" ./bin/frontend -prog -1 -role 2 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M "${m}" -OPT_BLOCK "${optB}" -repeats "${repeat}" -logFolder "${logFolder}"" &
wait;