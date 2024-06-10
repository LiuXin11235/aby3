"""Run basic performance benchmark tests.
"""

import time
import argparse
import json
import os
from utils import get_k

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"

# gtype_list = ["random", "star", "powerlaw", "bipartite", "tree"]
# gtype_list = ["random"]
gtype_list = ["random", "geometric", "powerlaw", "bipartite"]

format_configs = {
    "privGraph": {
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/privGraph/",
        "record_pattern": "(.*?)_n-(\d+)_k-(\d+)-(\d+)",
        "n": [1024],
        "e": -1,
        "k": [32],
        "n_stash_size": [32],
        "n_pack_size": [16],
        "e_stash_size": [1024],
        "e_pack_size": [32],
        "config_keys" : ["gtype", "n", "e", "k", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery"],
    },
    "adjmat":{
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/adjmat/",
        "record_pattern": "(.*?)_n-(\d+)-(\d+)",
        "n": [1024],
        "e": -1,
        "n_stash_size": [32],
        "n_pack_size": [16],
        "e_stash_size": [1024],
        "e_pack_size": [32],
        "config_keys" : ["gtype", "n", "e", "se", "pe", "sn", "pn"],
        "performance_keys": ["GraphLoad", "EdgeOramInit", "NodeOramInit", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery"],
    },
    "edgelist":{
        "prefix": MAIN_FOLDER + "/data/baseline/",
        "record_folder": MAIN_FOLDER + "/record/edgelist/",
        "record_pattern": "(.*?)_n-(\d+)-(\d+)",
        "n": [16384],
        "e": -1,
        "config_keys" : ["gtype", "n", "e"],
        "performance_keys": ["GraphLoad", "EdgeExistQuery", "OuttingEdgesCountQuery", "NeighborsGetQuery"],
    },
}

REPEAT_TIMES = 0

if __name__ == "__main__":
    
    # get the test configs. 
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="privGraph", help="target graph format")
    args = parser.parse_args()
    target = args.target
    
    # prepare the file. 
    os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")
    
    # run the privGraph profilings.
    target_config = format_configs[target] 
    
    if(not os.path.exists(target_config["record_folder"])):
        os.makedirs(target_config["record_folder"])
    if(not os.path.exists(target_config["prefix"])):
        os.makedirs(target_config["prefix"])
    e = -1
    
    if target == "privGraph":
    
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n, k = target_config["n"][i], target_config["k"][i]
                n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]

                # self define k.
                edge_list_meta_file = f"{target_config['prefix']}{gtype}_n-{n}_edge_list_meta.txt"
                if not os.path.exists(edge_list_meta_file):
                    print(f"File {edge_list_meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {target_config['prefix']}{gtype}_n-{n} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)

                with open(edge_list_meta_file, "r") as f:
                    edge_list_info = f.readline()
                    numbers = edge_list_info.split(" ")
                    V, E = int(numbers[0]), int(numbers[1])
                    assert V == n
                    k = get_k(V, E)
                
                # corresponsing data files.
                data_prefix = f"{gtype}_n-{n}_k-{k}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_meta.txt"
                graph_data_file = f"{file_prefix}_2dpartition.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --k {k}"
                    os.system(generate_command)
                
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -privGraph -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"
                    
                    print(f"Run the privGraph profiling with {run_args}")
                    os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
        
    # run the adjmat profilings.
    if target == "adjmat":
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n = target_config["n"][i]
                n_stash_size, n_pack_size = target_config["n_stash_size"][i], target_config["n_pack_size"][i]
                e_stash_size, e_pack_size = target_config["e_stash_size"][i], target_config["e_pack_size"][i]
                
                data_prefix = f"{gtype}_n-{n}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_edge_list_meta.txt"
                graph_data_file = f"{file_prefix}_edge_list.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)
                
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -adjmat true -prefix {data_prefix} -rcounter {count} -noram_stash_size {n_stash_size} -noram_pack_size {n_pack_size} -eoram_stash_size {e_stash_size} -eoram_pack_size {e_pack_size}"
                    
                    print(f"Run the adjmat profiling with {run_args}")
                    os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
        
    # run the edgelist profilings.
    if target == "edgelist":
        for gtype in gtype_list:
            for i in range(len(target_config["n"])):
                n = target_config["n"][i]
                
                data_prefix = f"{gtype}_n-{n}"
                file_prefix = f"{target_config['prefix']}{data_prefix}"
                meta_file = f"{file_prefix}_edge_list_meta.txt"
                graph_data_file = f"{file_prefix}_edge_list.txt"
                
                # if data do not exist, generate the data.
                if not os.path.exists(meta_file):
                    print(f"File {meta_file} does not exist. Generate the data...")
                    generate_command = f"python ./aby3-GraphQuery/privGraphQuery/micro_benchmark_generation.py --file_prefix {file_prefix} --type {gtype} --n {n} --e {e} --saving_type edgelist"
                    os.system(generate_command)
                
                # run the profiling.
                for j in range(REPEAT_TIMES):
                    count = j+1
                    run_args = f" -edgelist true -prefix {data_prefix} -rcounter {count}"
                    
                    print(f"Run the edgelist profiling with {run_args}")
                    os.system(f"./Eval/dis_exec.sh \"{run_args}\"")

                    # show some debug info. 
                    os.system(f"cat ./debug.txt; rm ./debug.txt")
                

    
    