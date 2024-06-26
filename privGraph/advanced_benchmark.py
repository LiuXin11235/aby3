"""Run advanced applications.
"""
import argparse
import json
import os
import random
from utils import get_k 

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"
record_folder = MAIN_FOLDER + "/record/adv_application/"
data_folder = ""
data_prefix = ""

n_stash_size, n_pack_size, e_stash_size, e_pack_size = 32, 16, 1024, 32

data_type = "real_world"

if(data_type == "synthetic"):
    data_folder = MAIN_FOLDER + "/data/baseline/"
    data_prefix = "bipartite_n-1024_k-4"
else:
    data_folder = MAIN_FOLDER + "/data/realworld/"
    data_prefix = "slashdot"
    
    
if __name__ == "__main__":
    
    # get the test configs.
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=str, default="cycle_detect", help="The target advanced applications")
    args = parser.parse_args()
    
    os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")
    
    if(not os.path.exists(data_folder)):
        os.makedirs(data_folder)
    if(not os.path.exists(record_folder)):
        os.makedirs(record_folder)
        
    if(args.target == "neighbor_stats"):
        # generate the property_graph.
        meta_file = data_folder + data_prefix + "_meta.txt"
        print(meta_file)
        with open(meta_file, "r") as f:
            line = f.readline()
            v, e, b, k, l = map(int, line.split())
        print(f"v: {v}, e: {e}, b: {b}, k: {k}, l: {l}")
        # generate the random property graph.
        property_file = data_folder + data_prefix + "_2dproperty.txt"
        with open(property_file, "w") as f:
            for i in range(b*b):
                for j in range(l):
                    random_number = random.randint(1, 10000)
                    f.write(f"{random_number}\n")
                
        
        edgelist_meta_file = data_folder + data_prefix + "_edge_list_meta.txt"
        with open(edgelist_meta_file, "r") as f:
            line = f.readline()
            v, e = map(int, line.split())
        print(f"v: {v}, e: {e}")
        
        # generate the random edgelist.
        property_edgelist_file = data_folder + data_prefix + "_edge_list_property.txt"
        with open(property_edgelist_file, "w") as f:
            for i in range(e):
                random_number = random.randint(1, 10000)
                f.write(f"{random_number}\n")
    
    # set the graph.
    run_args = f" -{args.target} -prefix {data_prefix} -rcounter 1 -noram_stash_size {n_stash_size} -eoram_stash_size {e_stash_size} -noram_pack_size {n_pack_size} -eoram_pack_size {e_pack_size} -data_folder {data_folder} -record_folder {record_folder}"
    
    print(run_args)
    print(f"Run {args.target}...")
    os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
    os.system(f"cat ./debug.txt; rm ./debug.txt")
    