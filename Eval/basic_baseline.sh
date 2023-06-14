task_list=("search" "select" "index" "average")
log_folder_list=(./Record/Record_search_base ./Record/Record_select_base ./Record/Record_index_base ./Record/Record_average_base)

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

# cleanup the old folder and create a clean one.
for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done

# compile
python build.py

# synchroonize with others
scp ./bin/frontend aby31:~/aby3/bin &
scp ./bin/frontend aby32:~/aby3/bin &
wait;

## Direct parallel baseline
# test settings.
N_list=(8388608 1048576)
M_list=(4 2 1 1 1 1)
simulate_task_num_list=(2 4 8 8 8 8)
optB_list=(262144)
repeat=3; test_times=5; retry_threshold=5
total_bw=10000; latency=0.003

# cleanup the old folder and create a clean one.
for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done

for (( i=0; i<${#task_list[@]}; i++ )); do

  # set the task the correspondng log folder
  task=${task_list[i]};
  log_folder=${log_folder_list[i]};

  for taskN in ${simulate_task_num_list[@]}; do

    # only set the bandwidth for simulation.
    test_bw=$(python -c "print($total_bw / ($taskN * 2) )")
    echo "test_bw = "${test_bw}
    # cleanup the settings on the eth
    ./Eval/network_clean.sh

    # setup the network
    ./Eval/network_set.sh ${test_bw} ${latency}

    for (( k=0; k<$test_times; k++ )); do
      for N in ${N_list[@]}; do
        for M in ${M_list[@]}; do
        # for M in ${M_list[@]}; do
          for optB in ${optB_list[@]}; do
            # run the tasks with retrying.
            j=0;
            while [ $j -lt $retry_threshold ]; do
              timeout 150m ./Eval/mpi_subtask.sh taskN=1 n=$N m=$M repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
              if [ $? -eq 0 ]; then
                break; 
              fi
              j=$(expr $j + 1);
              if [ $j -eq $retry_threshold ]; then
                echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
              fi
            done;
          done;
        done;
      done;
    done;
  done;

  # tar the results.
  tar cvf $log_folder.base.tar $log_folder;
done;