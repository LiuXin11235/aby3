prog=$1;
testFlag=$2;
if [ ! -n "${testFlag}"]; then
    testFlag=-1;
fi
echo "prog: "${prog}" | testFlag: "${testFlag};

if [ ${prog} == 0 ]; then
    keywords="basicOps";
elif [ ${prog} == 1 ]; then
    keywords="disBasicOps";
elif [ ${prog} == 2 ]; then
    keywords="test";
elif [ ${prog} == 3 ]; then
    keywords="cipher-index";
else
    keywords="invalid";
    echo "invalid prog";
fi

day=$(date +%m-%d);
timeStamp=$(date +"%H%M");

echo ${day}
echo ${timeStamp}


# compile
python build.py;

# synchronize with others
scp ./bin/frontend aby3-1:~/aby3/bin/ &
scp ./bin/frontend aby3-2:~/aby3/bin/ &
wait;

# set the log file
logFolder=./Record/Record${day}/
logFile=${logFolder}log-${keywords}-${day}${timeStamp}

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder};
fi

# run functions
./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} >> ${logFile} &
ssh aby3-1 "cd ./aby3/; ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" >> ./log" &
ssh aby3-2 "cd ./aby3/; ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" >> ./log" &
wait;