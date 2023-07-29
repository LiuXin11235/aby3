task="large_index"
log_folder=./Record/Record_index
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

N=10000000; M=1; repeat=1; optB=100000; taskN=16;
./Eval/mpi_subtask.sh taskN=$taskN n=$N m=$M repeat=$repeat task=$task optB=$optB log_folder=$log_folder/