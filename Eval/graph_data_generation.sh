# # generate test data.
# python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/micro_benchmark/tmp_graph;

# # generate one star graph for debugging.
# python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/micro_benchmark/star --type "star";

# generate the test adjancency matrix data.
# python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/micro_benchmark/adj_tmp --saving_type "edgelist" --type "random";

# generate the test adjancency matrix and edge list data.
graph_type_list=("random" "star" "powerlaw" "bipartite" "complete")
n_exp_list=(10 15)
for n_exp in ${n_exp_list[@]}; do
    n=`echo "2^$n_exp" | bc`;
    for gtype in ${graph_type_list[@]}; do
        python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix /root/aby3/aby3-GraphQuery/data/baseline/${gtype}_n=${n} --saving_type "edgelist" --type ${gtype} --n ${n} --e -1;
    done;
done;