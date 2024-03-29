# prepare the main
cp ./frontend/main.soram ./frontend/main.cpp

# log folder
logFolder=/root/aby3/Record/sqrt_oram/

if [ ! -d ${logFolder} ]; then
    mkdir ${logFolder}
fi

# compile the project
python build.py

# target tests.
M_list=(1048576 16777216)
N_list=(1 1)

# run the tests
task_args=" -BenchmarkORAM"

for((i=0; i<${#M_list[@]}; i++)); do
    M=${M_list[i]}; N=${N_list[i]};
    parames_args=" -M ${M} -N ${N}"
    logFile=${logFolder}log-config-M=${M}-N=${N}.log
    logging_folder_args=" -logFile ${logFile}"
    test_args=${task_args}${parames_args}${logging_folder_args}
    # echo ${test_args}
    ./Eval/basic/dis_exec.sh "${test_args}"
done