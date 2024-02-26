# generate test data.
python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/privGraphQuery/data/micro_benchmark/tmp_graph;

python build.py;

test_args=" -GraphQuery"
./Eval/dis_exec.sh "${test_args}"
wait;

cat ./debug.txt
rm ./debug.txt