args_list=$1
echo ${args_list}

# single-machine version
# ./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list} &
# ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list} &
# wait;


# # THE FOLLOWING IS FOR DISTRIBUTED TEST
# scp ./out/build/linux/frontend/frontend aby31:~/aby3/out/build/linux/frontend/ &
# scp ./out/build/linux/frontend/frontend aby32:~/aby3/out/build/linux/frontend/ &
# wait;

# distributed version
./out/build/linux/frontend/frontend -prog -1 -role 0 ${args_list} &
# public-private key authentication
# ssh aby31 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list}" &
# ssh aby32 "cd ./aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list}" &
# password authentication
sshpass -p aby3 ssh aby31 "cd ~/FedDetect/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 1 ${args_list}" &
sshpass -p aby3 ssh aby32 "cd ~/FedDetect/aby3/; ./out/build/linux/frontend/frontend -prog -1 -role 2 ${args_list}" &
wait;
