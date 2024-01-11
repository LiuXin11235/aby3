args_list=$1
echo ${args_list}
./bin/frontend -prog -1 -role 0 ${args_list} &
# ./bin/frontend -prog -1 -role 1 ${args_list} &
# ./bin/frontend -prog -1 -role 2 ${args_list} &
ssh aby31 "cd ./aby3/; ./bin/frontend -prog -1 -role 1 ${args_list}" &
ssh aby32 "cd ./aby3/; ./bin/frontend -prog -1 -role 2 ${args_list}" &
wait;
