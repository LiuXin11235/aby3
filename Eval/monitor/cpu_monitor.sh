report=$1; process_keyword=$2

PIDS=$(pgrep -f ${process_keyword} | tr '\n' ',' | sed 's/,$//')

watch -n 0.5 -t "ps -p ${PIDS} -o pid,ppid,cmd,%cpu,%mem >> ${report}" >> /dev/null