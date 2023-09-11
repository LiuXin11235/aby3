task_list=("index" "bio_metric" "average" "search" "select" "mean_distance")
log_folder_list=(./Record/Record_index_base ./Record/Record_bio_metric_base ./Record/Record_average_base ./Record/Record_search_base ./Record/Record_select_base ./Record/Record_mean_distance_base)


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
N_list=(1048576 4194304)
M_list=(8 4)
simulate_task_num_list=(2 4 8 16 16 16)
optB_list=(1048576 1048576 1048576 1048576 1048576 1048576)
repeat=1; test_times=1; retry_threshold=5
total_bw=20000; latency=0.003

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

  for (( t=0; t<${#simulate_task_num_list[@]}; t++ )); do
  # for taskN in ${simulate_task_num_list[@]}; do
    taskN=${simulate_task_num_list[t]}
    M=${M_list[i]}
    K=10

    # only set the bandwidth for simulation.
    test_bw=$(python -c "print($total_bw / ($taskN) )")
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
              timeout 15m ./Eval/mpi_subtask.sh taskN=1 n=$N m=$M k=$K repeat=$repeat task=$task optB=$optB log_folder=$log_folder/
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