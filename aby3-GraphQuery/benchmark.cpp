#include "benchmark.h"

using namespace oc;
using namespace aby3;

int privGraph_performance_profiling(oc::CLP& cmd){

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
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/profiling.igraph/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/privGraph/";
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

    // get the oram configurations, todo - define the parameters.
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
    if(role == 0) debug_info("Eoram construction...");
    timer.start("EdgeOramInit");
    secGraphEngine.edge_block_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    timer.start("NodeOramInit");
    secGraphEngine.node_edges_oram_initialization(noram_stash_size, noram_pack_size); 
    timer.end("NodeOramInit");

    if(role == 0) debug_info("Noram init success");

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
    aby3::si64Matrix out_edges = outting_edge_count(snode, snode_log_idx, secGraphEngine);
    timer.end("OuttingEdgesCountQuery");

    if(role == 0) debug_info("Outting edges count query success");

    // 3) neighbors get query.
    timer.start("NeighborsGetQuery");
    aby3::sbMatrix neighbors = outting_neighbors(snode, snode_log_idx, secGraphEngine);
    timer.end("NeighborsGetQuery");

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        secGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
    }

    return 0;
}

int adj_performance_profiling(oc::CLP& cmd){

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
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/baseline/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/adjmat/";
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

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the Mat Query Engine...");

    // get the oram configurations, todo - define the parameters.
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

    // get the timer.
    Timer& timer = Timer::getInstance();
    timer.clear_records();

    // graph loading.
    timer.start("GraphLoad");
    AdjGraphQueryEngine adjGraphEngine(party_info, meta_file, graph_data_file);
    timer.end("GraphLoad");

    if(role == 0) debug_info("Graph loaded successfully!");

    // oram construction.
    if(role == 0) debug_info("Eoram construction...");
    timer.start("EdgeOramInit");
    adjGraphEngine.edge_oram_initialization(eoram_stash_size, eoram_pack_size);
    timer.end("EdgeOramInit");

    if(role == 0) debug_info("Eoram init success\nNoram construction...");

    timer.start("NodeOramInit");
    adjGraphEngine.node_oram_initialization(noram_stash_size, noram_pack_size);
    timer.end("NodeOramInit");

    if(role == 0) debug_info("Noram init success");

    // 1) edge existence query.
    size_t v = adjGraphEngine.graph->v;
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(ps, role);
    boolIndex enode = boolIndex(pe, role);

    timer.start("EdgeExistQuery");
    boolShare flag = edge_existance(snode, enode, adjGraphEngine);
    timer.end("EdgeExistQuery");

    if(role == 0) debug_info("Edge existence query success");

    // 2) outting edges count query.
    timer.start("OuttingEdgesCountQuery");
    aby3::sbMatrix out_edges = outting_edge_count(snode, adjGraphEngine);
    timer.end("OuttingEdgesCountQuery");

    // 3) neighbors get query.
    timer.start("NeighborsGetQuery");
    aby3::sbMatrix neighbors = outting_neighbors(snode, adjGraphEngine);
    timer.end("NeighborsGetQuery");

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        adjGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
    }

    return 0;
}

int list_performance_profiling(oc::CLP& cmd){

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
    std::string graph_data_folder = "/root/aby3/aby3-GraphQuery/data/baseline/";
    std::string file_prefix = "tmp";
    std::string record_folder = "/root/aby3/aby3-GraphQuery/record/edgelist/";
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

    std::string meta_file = graph_data_folder + file_prefix + "_edge_list_meta.txt";
    std::string graph_data_file = graph_data_folder + file_prefix + "_edge_list.txt";
    std::string record_file = record_folder + file_prefix + "-" + std::to_string(record_counter) + ".txt";

    if (role == 0) debug_info("Profiling the Mat Query Engine...");

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    aby3Info party_info(role, enc, eval, runtime);

    if(role == 0) debug_info("Environment setup success");

    // get the timer.
    Timer& timer = Timer::getInstance();
    timer.clear_records();

    // graph loading.
    timer.start("GraphLoad");
    plainGraphList plainGraph(meta_file, graph_data_file);
    ListGraphQueryEngine listGraphEngine(party_info, plainGraph);
    timer.end("GraphLoad");

    if(role == 0) debug_info("Graph loaded successfully!");

    // 1) edge existence query.
    size_t v = plainGraph.v;
    size_t ps = 0, pe = v-1;
    boolIndex snode = boolIndex(ps, role);
    boolIndex enode = boolIndex(pe, role);
    timer.start("EdgeExistQuery");
    boolShare flag = edge_existance(snode, enode, listGraphEngine);
    timer.end("EdgeExistQuery");

    if(role == 0) debug_info("Edge existence query success");

    // 2) outting edges count query.
    timer.start("OuttingEdgesCountQuery");
    aby3::sbMatrix out_edges = outting_edge_count(snode, listGraphEngine);
    timer.end("OuttingEdgesCountQuery");

    // 3) neighbors get query.
    timer.start("NeighborsGetQuery");
    aby3::sbMatrix neighbors = outting_neighbors(snode, listGraphEngine);
    timer.end("NeighborsGetQuery");

    // print the timer records.
    if(role == 0){
        std::ofstream stream(record_file, std::ios::app);
        listGraphEngine.print_configs(stream);
        timer.print_total("milliseconds", stream);
    }


    return 0;
}