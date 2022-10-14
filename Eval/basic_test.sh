# compile
python build.py;

# synchronize with others
scp ./bin/frontend aby3-130:~/aby3/bin/ &
scp ./bin/frontend aby3-131:~/aby3/bin/ &
wait;

# set the log file
logFolder=./Record0927/
logFile=${logFolder}log-basicOps

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder};
fi

rm -r ${logFolder}*

# run functions
./bin/frontend -role 0 >> ${logFile} &
ssh aby3-130 "cd ./aby3/; ./bin/frontend -role 1 >> ./log" &
ssh aby3-131 "cd ./aby3/; ./bin/frontend -role 2 >> ./log" &
wait;