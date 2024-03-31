report=$1; process_keyword=$2

sudo nethogs -d 0.5 -t | grep ${process_keyword} >> ${report}