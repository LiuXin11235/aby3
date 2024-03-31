# prepare the main file.
cp ./frontend/main.pta ./frontend/main.cpp
USER_FOLDER=/home/tsingj_ubuntu/fanxy/PtA/aby3

# generate the profiling metrics
python $USER_FOLDER/PtA_deploy/generate_execution_shells.py

# compile the project
python build.py

# run the generated shell scripts
${USER_FOLDER}/Eval/pta_execution.sh