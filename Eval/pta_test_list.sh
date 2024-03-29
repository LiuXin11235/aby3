task_list=("cipher_index")

logFolder=/root/aby3/Record/
ptaLogFolder=${logFolder}pta/

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder}
fi

if [ ! -d ${ptaLogFolder} ]; then
    mkdir ${ptaLogFolder}
fi

declare -A M_lists
M_lists["cipher_index"]="4194304 4194304 4194304"

declare -A N_lists
N_lists["cipher_index"]="1 1 1"

declare -A c_lists
c_lists["cipher_index"]="5 10 20"

declare -A optB_lists
optB_lists["cipher_index"]="32768 65536 131072"

declare -A logFolders
logFolders["cipher_index"]=${ptaLogFolder}"cipher_index/"


# compile the main.
cp ./frontend/main.pta ./frontend/main.cpp
python build.py

for task in ${task_list[@]}; do

    if [ ! -d ${logFolders[$task]} ]; then
        mkdir ${logFolders[$task]}
    fi

    # get the parameters
    read -a M_list <<< ${M_lists[$task]}; 
    read -a N_list <<< ${N_lists[$task]};
    read -a c_list <<< ${c_lists[$task]};
    read -a optB_list <<< ${optB_lists[$task]};

    param_len=${#M_list[@]}

    echo "param_len="$param_len

    task_args=" -${task}"

    for ((i=0; i<param_len; i++)); do
        N=${N_list[i]}; M=${M_list[i]}; c=${c_list[i]}; optB=${optB_list[i]};
        parames_args=" -N ${N} -M ${M} -c ${c} -B ${optB}"
        # get the logging file.
        logFile=${logFolders[$task]}log-config-N=${N}-M=${M}-c=${c}-B=${optB}.log
        logging_folder_args=" -logFile ${logFile}"

        test_args=${task_args}${parames_args}${logging_folder_args} 

        ./Eval/basic/dis_mpi_exec.sh ${c} "${test_args}"

    done;
done