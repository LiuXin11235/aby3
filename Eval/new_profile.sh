task="profile_index"

total_bw=10000
available_cores_list=(256 128 64 32 16 8 4 2 1)
latency=0.03

# compile.
python build.py

# synchronize with others.
scp ./bin/frontend aby31:~/aby3/bin &
scp ./bin/frontend aby32:~/aby3/bin &
wait;

for available_cores in ${available_cores_list[@]}; do

  test_bw=$(python -c "print($total_bw / ($available_cores * 2) )")
  echo "test_bw = "${test_bw}

  # cleanup the settings on the eth
  ./Eval/network_clean.sh

  # add a new, clean log folder.
  log_folder=/root/aby3/Record/Record_probe

  day=$(date +%m-%d);
  timeStamp=$(date +"%H%M%s");

  # backup the last log folder.
  if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
  else
    rm -r ${log_folder};
    mkdir ${log_folder};
  fi

  # setup the network
  ./Eval/network_set.sh ${test_bw} ${latency}

  # evaluate on aby3
  n=1; optB=64; m=1000000000; repeat=100; epsilon=1000; gap=100

  timeout 100m ./Eval/mpi_subtask.sh taskN=1 n=${n} m=${m} repeat=${repeat} task=${task} optB=${optB} log_folder=${log_folder}/ epsilon=${epsilon} gap=${gap} vec_start=${optB}

  cp ./Record/Record_probe/probe.log ./Record/probe_info/probe_${test_bw}bw.log
  cp ./Record/Record_probe/probe.res ./Record/probe_info/probe_${test_bw}bw.res

done;
