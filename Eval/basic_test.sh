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
elif [ ${prog} == 4 ]; then
    keywords="new-api";
elif [ ${prog} == 5 ]; then
    keywords="ptr";
else
    keywords="invalid";
    echo "invalid prog";
fi

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

echo ${day}
echo ${timeStamp}

rm /root/aby3/debug.txt

# compile
python build.py;

# synchronize with others
scp ./bin/frontend aby31:~/aby3/bin/ &
scp ./bin/frontend aby32:~/aby3/bin/ &
wait;

# set the log file
logFolder=./Record/Record${day}/
logFile=${logFolder}log-${keywords}-${day}${timeStamp}

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder};
fi

# run functions
./bin/frontend -prog ${prog} -role 0 -testFlag ${testFlag} >> ${logFile} &
ssh aby31 "cd ./aby3/; ./bin/frontend -prog "${prog}" -role 1 -testFlag "${testFlag}" >> ./log" &
ssh aby32 "cd ./aby3/; ./bin/frontend -prog "${prog}" -role 2 -testFlag "${testFlag}" >> ./log" &
wait;

cat $logFile