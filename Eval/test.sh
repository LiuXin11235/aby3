# prepare the test cpp.
cp ./frontend/main.test ./frontend/main.cpp

# compile the main.
python build.py

# clean debugging files party-*.txt if exist.
for pfile in ./party-*.txt; do
    rm ${pfile};
done

# # synchronize with others
# scp ./bin/frontend aby31:~/aby3/bin
# scp ./bin/frontend aby32:~/aby3/bin
# wait;

# run the tests
# current tests: 
# 1) -Bool : boolean share tests; 
# 2) -Arith : arithmetic share tests; 
# 3) -ORAM : ORAM tests; 
# 4) -Init : initialization tests, including the correlated shares; 
# 5) -Shuffle : secure shuffling tests.
test_args=" -ORAM"
./Eval/dis_exec.sh "${test_args}"
wait;

# scp aby31:~/aby3/debug.txt ./debug-p1.txt
# scp aby32:~/aby3/debug.txt ./debug-p2.txt

cat ./debug.txt
rm ./debug.txt
