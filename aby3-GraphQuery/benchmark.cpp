#include "benchmark.h"

using namespace oc;
using namespace aby3;

int performance_profiling(oc::CLP& cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    // get the graph file parameters.
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/profiling/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/";
    int record_counter = -1;

    if(cmd.isSet("prefix")){
        auto keys = cmd.getMany<std::string>("prefix");
        file_prefix = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("prefix must be set!");
    }

    if(cmd.isSet("rcounter")){
        auto keys = cmd.getMany<int>("rcounter");
        record_counter = keys[0];
    }
    else{
        THROW_RUNTIME_ERROR("rcounter must be set!");
    }


    std::string meta_file = graph_data_folder + file_prefix + "_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_2dpartition.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the privGraph Query Engine...");

    // get the oram configurations.
    size_t noram_stash_size = 1 << 4;
    size_t noram_pack_size = 1 << 2;
    size_t eoram_stash_size = 1 << 8;
    size_t eoram_pack_size = 1 << 4;

    if(cmd.isSet("noram_stash_size")){
        auto keys = cmd.getMany<size_t>("noram_stash_size");
        noram_stash_size = keys[0];
    }
    else{
        debug_info("noram_stash_size is not set, using default value: " + std::to_string(noram_stash_size));
    }
    if(cmd.isSet("noram_pack_size")){
        auto keys = cmd.getMany<size_t>("noram_pack_size");
        noram_pack_size = keys[0];
    }
    else{
        debug_info("noram_pack_size is not set, using default value: " + std::to_string(noram_pack_size));
    }
    if(cmd.isSet("eoram_stash_size")){
        auto keys = cmd.getMany<size_t>("eoram_stash_size");
        eoram_stash_size = keys[0];
    }
    else{
        debug_info("eoram_stash_size is not set, using default value: " + std::to_string(eoram_stash_size));
    }
    if(cmd.isSet("eoram_pack_size")){
        auto keys = cmd.getMany<size_t>("eoram_pack_size");
        eoram_pack_size = keys[0];
    }
    else{
        debug_info("eoram_pack_size is not set, using default value: " + std::to_string(eoram_pack_size));
    }

    if(role == 0) debug_info("Getting oram configs success");

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    if(role == 0) debug_info("Environment setup success");

    // load the graph.
    Timer& timer = Timer::getInstance();
    timer.start("GraphLoad");
    GraphQueryEngine secGraphEngine(party_info, meta_file, graph_data_file);
    timer.end("GraphLoad");

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    timer.start("EdgeOramInit");
    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");

    if(role == 0) debug_info("eoram init success");

    timer.start("NodeOramInit");
    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size); 
    timer.end("NodeOramInit");

    if(role == 0) debug_info("noram init success");

    if(role == 0){
        std::ofstream ofs(debugFile, std::ios_base::app);
        secGraphEngine.print_configs(ofs);
        ofs.close();
    }

    size_t b = secGraphEngine.graph->b;
    size_t b2 = secGraphEngine.graph->edge_list_size;
    size_t v = secGraphEngine.graph->v;

    // edge block fetch
    for(int i=0; i<eoram_stash_size; i++){
        boolIndex tar_ind = boolIndex((i % b2), role);
        std::string timer_key = timer.get_key("EdgeBlockFetch");
        timer.start(timer_key);
        secGraphEngine.get_edge_block(tar_ind);
        timer.end(timer_key);
    }

    if(role == 0) debug_info("Edge block fetch success");

    // node edges block fetch
    for(int i=0; i<noram_stash_size; i++){
        boolIndex tar_ind = boolIndex((i % b), role);
        std::string timer_key = timer.get_key("NodeEdgesBlockFetch");
        timer.start(timer_key);
        secGraphEngine.get_node_edges(tar_ind);
        timer.end(timer_key);
    }

    if(role == 0) debug_info("node block fetch success");
    
    // reinit secGraphEngine.
    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);
    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size);

    if(role == 0) debug_info("reinit success");

    // basic graph query process.
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(secGraphEngine.get_block_index(ps), role);
    boolIndex enode = boolIndex(secGraphEngine.get_block_index(pe), role);
    boolIndex edge_log_idx = boolIndex(secGraphEngine.get_edge_block_index(ps, pe), role);
    boolIndex snode_log_idx = boolIndex(secGraphEngine.get_block_index(ps), role);

    // 1) edge existence query.
    timer.start("EdgeExistQuery");
    boolShare flag = edge_existance(snode, enode, edge_log_idx, secGraphEngine);
    timer.end("EdgeExistQuery");

    if(role == 0) debug_info("Edge existence query success");

    // 2) outting edges count query.
    timer.start("OuttingEdgesCountQuery");
    aby3::sbMatrix out_edges = outting_edge_count(snode, snode_log_idx, secGraphEngine);
    timer.end("OuttingEdgesCountQuery");

    if(role == 0) debug_info("Outting edges count query success");

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
        // timer.print_records("EdgeBlockFetch", "milliseconds", stream);
        // timer.print_records("NodeEdgesBlockFetch", "milliseconds", stream);
    }

    return 0;
}