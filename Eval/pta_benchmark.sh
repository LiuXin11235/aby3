# prepare the main file.
cp ./frontend/main.pta ./frontend/main.cpp

# generate the profiling metrics
python /root/aby3/PtA_deploy/generate_execution_shells.py

# compile the project
python build.py

# run the generated shell scripts
/root/aby3/Eval/pta_execution.sh