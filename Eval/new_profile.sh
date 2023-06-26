task_list=("prof_bio_metric" "prof_average" "prof_mean_distance" "prof_cipher_index" "prof_sort" "prof_rank" "prof_max")
log_folder_list=(./Record/Prof_bio_metric ./Record/Prof_average ./Record/Prof_mean_distance ./Record/Prof_index ./Record/Prof_sort ./Record/Prof_rank ./Record/Prof_max)

# task_list=("prof_bio_metric" "prof_max" "prof_sort" "prof_rank")
# log_folder_list=(./Record/Prof_bio_metric ./Record/Prof_max ./Record/Prof_sort ./Record/Prof_rank)

# cleanup the old folder and create a clean one.
for log_folder in ${log_folder_list[@]}; do
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi
done

# construct the res folder
res_folder=./Record/probe_info
if [ ! -d ${res_folder} ]; then
  mkdir ${res_folder}
else
  rm -r ${res_folder};
  mkdir ${res_folder};
fi


total_bw=10000
# available_cores_list=(256 64 16 4 1)
available_cores_list=(256)
latency=0.03

# compile.
python build.py

# # synchronize with others.
# scp ./bin/frontend aby31:~/aby3/bin &
# scp ./bin/frontend aby32:~/aby3/bin &
# wait;

retry_threshold=5;

for (( i=0; i<${#task_list[@]}; i++ )); do
  task=${task_list[i]}; log_folder=${log_folder_list[i]};

  for available_cores in ${available_cores_list[@]}; do

    test_bw=$(python -c "print($total_bw / ($available_cores * 2) )")
    echo "test_bw = "${test_bw}

    # cleanup the settings on the eth
    # ./Eval/network_clean.sh

    # add a new, clean log folder.
    day=$(date +%m-%d);
    timeStamp=$(date +"%H%M%s");

    # setup the network
    # ./Eval/network_set.sh ${test_bw} ${latency}

    # evaluate on aby3
    n=1; optB=1024; m=1000000000; repeat=10; epsilon=0.1; gap=100; K=2;

    j=0;
    while [ $j -lt $retry_threshold ]; do
      timeout 10m ./Eval/eval_subtask.sh taskN=1 n=${n} m=${m} repeat=${repeat} task=${task} optB=${optB} log_folder=${log_folder}/ epsilon=${epsilon} gap=${gap} vec_start=${optB} k=${K};
      if [ $? -eq 0 ]; then
        break; 
      fi
      j=$(expr $j + 1);
      if [ $j -eq $retry_threshold ]; then
        echo "Max retry: "${n}-${taskN} >> ${log_folder}/error.log;
      fi
    done;

    cp ${log_folder}/probe.log ${res_folder}/probe-task_${task}_${test_bw}bw.log
    cp ${log_folder}/probe.res ${res_folder}/probe-task_${task}_${test_bw}bw.res
  done;
done;
