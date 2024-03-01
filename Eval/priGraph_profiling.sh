# define the main file.
cp ./frontend/main.pgp ./frontend/main.cpp

# build the frontend.exe
# python build.py

# config_params.
# gtype_list=("random" "complete" "star" "powerlaw" "bipartite")
gtype_list=("random" "star" "powerlaw" "bipartite")

for gtype in ${gtype_list[@]}; do
    nExp=10; n=`echo "2^$nExp" | bc`;
    kExp=5; k=`echo "2^$kExp" | bc`;
    e=-1

    folder_prefix="/root/aby3/aby3-GraphQuery/data/profiling/"
    prefix="${gtype}_n-${n}_k-${k}"
    file_prefix="${folder_prefix}${prefix}"

    # generate data.
    # check whether data exist.
    meta_file=${file_prefix}_meta.txt
    graph_data_file=${file_prefix}_2dpartition.txt

    if [ -f ${meta_file} -a -f ${graph_data_file} ]; then
        echo "Data ${prefix} already exists. Skip the data generation."
    else
        echo "Data ${prefix} does not exist. Generate the data..."
        python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix ${file_prefix} --type ${gtype} --n ${n} --e ${e} --k ${k}
    fi

    echo "Data ${prefix} generation DONE."

done;

# n_stashExp=5; n_packExp=4;
# e_stashExp=10; e_packExp=5;
# n_stash=`echo "2^$n_stashExp" | bc`; n_pack=`echo "2^$n_packExp" | bc`;
# e_stash=`echo "2^$e_stashExp" | bc`; e_pack=`echo "2^$e_packExp" | bc`;

# # run the profiling test.
# run_args=" -prefix ${prefix} -rcounter 0 -noram_stash_size ${n_stash} -noram_pack_size ${n_pack} -eoram_stash_size ${e_stash} -eoram_pack_size ${e_pack}"
# echo ${run_args}

# ./Eval/dis_exec.sh "${run_args}"

# cat ./debug.txt
# rm ./debug.txt