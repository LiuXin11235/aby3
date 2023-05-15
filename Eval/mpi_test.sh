prog=$1;
TASKS=20;
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

# run functions
mpirun -n ${TASKS} ./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} >> ${logFile} &
ssh aby31 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" >> ./log" &
ssh aby32 "cd ./aby3/; mpirun -n "${TASKS}" ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" >> ./log" &
wait;

cat $logFile;