task="profile_index"

total_bw=10000
available_cores=1024
test_bw=$(python -c "print($total_bw / ($available_cores * 2) )")
echo "test_bw = "${test_bw}
latency=0.03

# compile.
python build.py

# synchronize with others.
scp ./bin/frontend aby31:~/aby3/bin/ &
scp ./bin/frontend aby32:~/aby3/bin/ &
wait;

# add a new, clean log folder.
log_folder=/root/aby3/Record/Record_probe

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# backup the last log folder.
if [ ! -d ${log_folder} ]; then
  mkdir ${log_folder}
else
  # cp -r ${log_folder} ${log_folder}${day}${timeStamp};
  rm -r ${log_folder};
  mkdir ${log_folder};
fi

# # cleanup the settings on the eth
tc qdisc del dev eth0 root; tc qdisc del dev eth0 ingress;
ssh aby31 "tc qdisc del dev eth0 root; tc qdisc del dev eth0 ingress";
ssh aby32 "tc qdisc del dev eth0 root; tc qdisc del dev eth0 ingress";

# setup the network
tc qdisc add dev eth0 handle 1: ingress;
tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1;
tc qdisc add dev eth0 root tbf rate ${test_bw}mbit latency ${latency}ms burst ${test_bw}mbit;
ssh aby31 "tc qdisc add dev eth0 handle 1: ingress; tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1; tc qdisc add dev eth0 root tbf rate "${test_bw}"mbit latency "${latency}"ms burst "${test_bw}"Mb";
ssh aby32 "tc qdisc add dev eth0 handle 1: ingress; tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1; tc qdisc add dev eth0 root tbf rate "${test_bw}"mbit latency "${latency}"ms burst "${test_bw}"Mb";

# evaluate on aby3
n=1; optB=32; m=1000000000; repeat=1000; epsilon=1000; gap=100

timeout 100m ./Eval/mpi_subtask.sh taskN=1 n=${n} m=${m} repeat=${repeat} task=${task} optB=${optB} log_folder=${log_folder}/ epsilon=${epsilon} gap=${gap} vec_start=${optB}

cp ./Record/Record_probe/probe.log ./Record/probe_info/probe_${test_bw}bw.log

# # # cleanup the settings on the eth
# tc qdisc del dev eth0 root
# ssh aby31 "tc qdisc del dev eth0 root";
# ssh aby32 "tc qdisc del dev eth0 root";


# # bandwidth_list=(1 5 100 1000)
# # bandwidth_list=(1)
# bandwidth_list=(20)
# latency=0.03

# # bw=100000
# # latency_list=(5 50 100)
# # latency_list=(1000 100 50 5 1)

# # for latency in ${latency_list[@]}; do
# for bw in ${bandwidth_list[@]}; do

#   # log_folder=/root/aby3/Record/Record_latency_${latency}Ms
#   log_folder=/root/aby3/Record/Record_bandwidth_${bw}Mb

#   day=$(date +%m-%d);
#   timeStamp=$(date +"%H%M%s");

#   # backup the last log folder.
#   if [ ! -d ${log_folder} ]; then
#     mkdir ${log_folder}
#   else
#     cp -r ${log_folder} ${log_folder}${day}${timeStamp};
#     rm -r ${log_folder};
#     mkdir ${log_folder};
#   fi

#   # compile
#   python build.py

#   # synchroonize with others
#   scp ./bin/frontend aby31:~/aby3/bin/ &
#   scp ./bin/frontend aby32:~/aby3/bin/ &
#   wait;

#   ./Eval/kill_all.sh frontend;

#   # tasks_list=(2 4 8 16 32 64)
#   tasks_list=(1)
#   n_list=(10 46 215 1000 4641 21544 100000 464158 2154434 10000000)
#   repeats_list=(1 1 1 1 1 1 1 1 1 1)
#   m=1

#   # # setup latency
#   # tc qdisc add dev eth0 root netem delay ${latency}ms;
#   # ssh aby31 "tc qdisc add dev eth0 root netem delay "${latency}"ms";
#   # ssh aby32 "tc qdisc add dev eth0 root netem delay "${latency}"ms";

#   # setup bandwidth
#   tc qdisc add dev eth0 handle 1: ingress;
#   tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1;
#   tc qdisc add dev eth0 root tbf rate ${bw}mbit latency ${latency}ms burst ${bw}mbit;
#   ssh aby31 "tc qdisc add dev eth0 handle 1: ingress; tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1; tc qdisc add dev eth0 root tbf rate "${bw}"mbit latency "${latency}"ms burst "${bw}"Mb";
#   ssh aby32 "tc qdisc add dev eth0 handle 1: ingress; tc filter add dev eth0 parent 1: protocol ip prio 50 u32 match ip src 0.0.0.0/0 flowid :1; tc qdisc add dev eth0 root tbf rate "${bw}"mbit latency "${latency}"ms burst "${bw}"Mb";


#   for taskN in ${tasks_list[@]}; do
#     for (( i=0; i<${#n_list[@]}; i++ )); do
#       optB=${n_list[i]};
#       # n=${optB}
#       n=$(($optB * $taskN))
#       repeat=${repeats_list[i]};

#       day=$(date +%m-%d);
#       timeStamp=$(date +"%H%M%s");
#       echo ${day}-${timeStamp}-${taskN}-${optB} >> ${log_folder}/log.time

#       j=0;
#       while [ $j -lt 5 ]; do
#         timeout 100m ./Eval/general_mpi.sh ${taskN} ${n} ${m} ${repeat} ${task} ${optB} ${log_folder}/
#         if [ $? -eq 0 ];then
#           break;
#         fi
#         echo ">>>>>>>> in retry!"${taskN}-${optB};
#         ./Eval/kill_all.sh frontend;
#         sleep 1m;
#         echo ">>>>>>>> after sleep!"${taskN}-${optB};
#         ipcs -m shmid;
#         echo "retry: "${j}" for "${optB}-${taskN} >> ${log_folder}/error.log;

#         rm ${log_folder}/log-config-N=${n}-TASKS=${taskN}-Vec=${optB}-0;
#         # for (( id=0; id<$taskN; id++ )); do
#         #   rm ${log_folder}/log-config-N=${n}-TASKS=${taskN}-Vec=${optB}-${id};
#         # done

#         j=$(expr $j + 1);
#         if [ $j -eq 5 ];then
#           echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
#         fi
#       done;

#       ./Eval/kill_all.sh frontend;
#       sleep 10;
#       ipcs -m shmid;
#     done;
#   done;

#   tar cvf ${log_folder}.tar ${log_folder};

#   # cleanup the settings on the eth
#   tc qdisc del dev eth0 root
#   ssh aby31 "tc qdisc del dev eth0 root";
#   ssh aby32 "tc qdisc del dev eth0 root";

# done;



# # for taskN in ${tasks_list[@]}; do
# #   for (( i=0; i<${#n_list[@]}; i++ )); do
# #     n=${n_list[i]};
# #     repeat=${repeats_list[i]};
    
# #     day=$(date +%m-%d);
# #     timeStamp=$(date +"%H%M%s");
# #     echo ${day}-${timeStamp}-${taskN}-${n} >> ${log_folder}/log.time

# #     mpirun -n ${taskN} ./bin/frontend -prog -1 -role 0 -N ${n} -TASK_NUM ${taskN} -FUNC ${task} -M 0 -OPT_BLOCK -1 -repeats ${repeat} &
# #     ssh aby31 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 1 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" >> ./log" &
# #     ssh aby32 "cd ./aby3/; mpirun -n "${taskN}" ./bin/frontend -prog -1 -role 2 -N "${n}" -TASK_NUM "${taskN}" -FUNC "${task}" -M 0 -OPT_BLOCK -1 -repeats "${repeat}" >> ./log" &
# #     wait;

# #     sleep 60;

# #   done;
# # done;

