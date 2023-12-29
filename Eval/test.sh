# compile the main.
python build.py

# rm ./debug.txt

# synchronize with others
scp ./bin/frontend aby31:~/aby3/bin
scp ./bin/frontend aby32:~/aby3/bin
wait;

test_args=" -Test true"
./Eval/dis_exec.sh "${test_args}"
wait;

cat ./debug.txt
rm ./debug.txt