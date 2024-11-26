# merge the data
cat lr_train/train_data_party*.csv > lr_train/train_data.csv
cat lr_train/train_label_party*.csv > lr_train/train_label.csv
cat lr_train/test_data_party*.csv > lr_train/test_data.csv
cat lr_train/test_label_party*.csv > lr_train/test_label.csv

# compile the main.
cp ./frontend/main.test ./frontend/main.cpp
current_path=$(pwd)
debugFile="${current_path}/debug.txt"
graphFolder="${current_path}/aby3-GORAM/data/"
echo "Current path: ${debugFile}"
python build.py --DEBUG_FILE ${debugFile} --GRAPH_FOLDER ${graphFolder}

# clean debugging files party-*.txt if exist.
for pfile in ./party-*.txt; do
    rm ${pfile};
done

# ./Eval/graph_test.sh

# # synchronize with others
# scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
# scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
# wait;

# run the tests
# current tests: 
# 1) -Bool : boolean share tests; 
# 2) -Arith : arithmetic share tests; 
# 3) -ORAM : ORAM tests; 
# 4) -Init : initialization tests, including the correlated shares; 
# 5) -Shuffle : secure shuffling tests.
# 6) -Graph : basic graph loading tests.
# 7) -GraphQuery : basic graph query tests (block fetching, edge exist & outting edges count.)
# 8) -Comm : test inter-party communication.
# 9) -Sort : test the sort functions.
# test_args=" -ORAM -Graph"
# test_args=" -Graph -GraphQuery -Sort"
# test_args=" -Sort -GraphQuery"
# test_args=" -GraphQuery"
# test_args=" -Bool -Comm -Graph -GraphQuery"
# test_args=" -Comm -Bool -Graph -GraphQuery"
# test_args=" -Shuffle -ORAM -Graph -GraphQuery -Sort"
# test_args=" -LR -Shuffle -Sort"
test_args=" -LR"
./Eval/dis_exec.sh "${test_args}"
wait;

cat ./debug.txt
rm ./debug.txt
