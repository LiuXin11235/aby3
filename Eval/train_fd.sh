# compile the main.
cp ./frontend/main.test ./frontend/main.cpp
current_path=$(pwd)
debugFile="${current_path}/debug.txt"
graphFolder="${current_path}/aby3-GORAM/data/"
echo "Current path: ${debugFile}"
python build.py --DEBUG_FILE ${debugFile} --GRAPH_FOLDER ${graphFolder}
cat lr_train_party/train_data* > lr_train/train_data.csv
cat lr_train_party/train_label* > lr_train/train_label.csv
cat lr_train_party/val_data* > lr_train/val_data.csv
cat lr_train_party/val_label* > lr_train/val_label.csv
cat lr_train_party/test_data* > lr_train/test_data.csv
cat lr_train_party/test_label* > lr_train/test_label.csv

# clean debugging files party-*.txt if exist.
for pfile in ./party-*.txt; do
    rm ${pfile};
done

# ./Eval/graph_test.sh

# synchronize with others
scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
wait;

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
test_args=" -LR_train"
./Eval/dis_exec.sh "${test_args}"
wait;

cat ./debug.txt
rm ./debug.txt

