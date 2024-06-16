import argparse
import pandas as pd
import numpy as np
import os

MAIN_FOLDER = "/root/aby3/aby3-GraphQuery"
real_world_data_folder = MAIN_FOLDER + "/data/realworld/"  
main_record_folder = MAIN_FOLDER + "/record/realworld/"

n_stash_size, n_pack_size, e_stash_size, e_pack_size = 32, 16, 1024, 32
test_format = ["privGraph", "adjmat", "edgelist"]

REPEAT_TIMES = 0

if __name__ == "__main__":

    # prepare the file.
    os.system("cp ./frontend/main.pgp ./frontend/main.cpp; python build.py")

    # get the test configs.
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="", help="target real-world graph")
    args = parser.parse_args()
    target = args.target

    for gformat in test_format:
        # different graph format.
        target_record_folder = main_record_folder + f"{gformat}/"
        if(not os.path.exists(target_record_folder)):
            os.makedirs(target_record_folder)

        run_args = f" -privGraph -prefix {target} -rcounter 1 -noram_stash_size {n_stash_size} -eoram_stash_size {e_stash_size} -noram_pack_size {n_pack_size} -eoram_pack_size {e_pack_size} -data_folder {real_world_data_folder} -record_folder {target_record_folder}"

        print(run_args)
        os.system(f"./Eval/dis_exec.sh \"{run_args}\"")
        os.system(f"cat ./debug.txt; rm ./debug.txt")
