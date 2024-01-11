# compile the main.
python build.py

# synchronize with others
scp ./bin/frontend aby31:~/aby3/bin
scp ./bin/frontend aby32:~/aby3/bin
wait;

test_args=" -Bool -Arith"
# test_args="-u " # origional aby3 test, failed everywhere.
./Eval/dis_exec.sh "${test_args}"
wait;

# scp aby31:~/aby3/debug.txt ./debug-p1.txt
# scp aby32:~/aby3/debug.txt ./debug-p2.txt

cat ./debug.txt
rm ./debug.txt