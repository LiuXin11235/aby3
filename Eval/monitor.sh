log_folder=./Record/Record_vector/log-monitor

day=$(date +%m-%d);
timeStamp=$(date +"%H%M%s");

while true; do
  echo ${day}-${timeStamp} >> ${log_folder}.cpu
  echo ${day}-${timeStamp} >> ${log_folder}.net
  top -b >> ${log_folder}.cpu 
  iftop -t >> ${log_folder}.net
done;