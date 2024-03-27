# generate test data.
python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/micro_benchmark/tmp_graph;

# generate one star graph for debugging.
python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/micro_benchmark/star --type "star";