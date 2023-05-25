prog=$1;
testFlag=-1

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# set the log file
logFolder=./MPIRecord/Record${day}/
logFile=${logFolder}log-${keywords}-${day}${timeStamp}

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder};
fi

# compile
python build.py;

# synchronize with others
scp ./bin/frontend aby31:~/aby3/bin/ &
scp ./bin/frontend aby32:~/aby3/bin/ &
wait;

N_LIST=(100000000)
M_LIST=(1)
TASKS_LIST=(128)
OPT_BLOCKS=(1000 10000 100000 500000 1000000 5000000)


for N in ${N_LIST[@]}; do
    for M in ${M_LIST[@]}; do
        for TASKS in ${TASKS_LIST[@]}; do
            for OPT_B in ${OPT_BLOCKS[@]}; do
                echo ${TASKS}
                # run functions
                mpirun -np 128 ./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} -N ${N} -M ${M} -TASK_NUM 128 -OPT_BLOCK ${OPT_B} >> ${logFile} &

                ssh aby31 "cd ./aby3/; mpirun -np 128 ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM 128 -OPT_BLOCK "${OPT_B}" >> ./log" &

                ssh aby32 "cd ./aby3/; mpirun -np 128 ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM 128 -OPT_BLOCK "${OPT_B}" >> ./log" &


                # mpirun -n ${TASKS} ./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} -N ${N} -M ${M} -TASK_NUM ${TASKS} -OPT_BLOCK ${OPT_B} >> ${logFile} &

                # ssh aby31 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" >> ./log" &

                # ssh aby32 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" -N "${N}" -M "${M}" -TASK_NUM "${TASKS}" -OPT_BLOCK "${OPT_B}" >> ./log" &
                wait;
            done;
        done;
    done;
done;
wait;

# run functions
# mpirun -n ${TASKS} ./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} >> ${logFile} &
# ssh aby31 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" >> ./log" &
# ssh aby32 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" >> ./log" &

# cat $logFile;