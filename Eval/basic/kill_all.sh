keyword=$1

pkill -f ${keyword}
ssh sosp1 "pkill -f ${keyword}"
ssh sosp2 "pkill -f ${keyword}"