#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include "pGraph.h"
#include "aby3-Basic/SqrtOram.h"
#include "aby3-Basic/Sort.h"

struct aby3Info{
    int pIdx;
    aby3::Sh3Encryptor *enc;
    aby3::Sh3Evaluator *eval;
    aby3::Sh3Runtime *runtime;

    aby3Info(int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime):
        pIdx(pIdx), enc(&enc), eval(&eval), runtime(&runtime){}
};

struct Graph2d {
    // each 2d-partation graph includes one node chuck list (length b) and one edge block list (length b^2).
    size_t v, e;
    size_t b, k, l;
    size_t edge_list_size;
    // the node_chuck_list is a list of sbMatrix, each sbMatrix is a list of k nodes.
    std::vector<aby3::sbMatrix> node_chuck_list;

    // the edge_block_list is a list of sbMatrix, each sbMatrix is a list of 2*l node tags (the first l is starting nodes and the last l is ending nodes).
    std::vector<aby3::sbMatrix> edge_block_list;

    // the node_edges_list is a list of sbMatrix, each sbMatrix is a list of 2*b*l node tags (the first b*l is starting nodes and the last b*l is ending nodes).
    std::vector<aby3::sbMatrix> node_edges_list;

    Graph2d(){}; // default constructor

    Graph2d(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info &party_info){
        plainGraph2d plain_graph(meta_data_file, edge_block_file);
        v = plain_graph.v;
        e = plain_graph.e;
        b = plain_graph.b;
        k = plain_graph.k;
        l = plain_graph.l;
        edge_list_size = plain_graph.edge_list_size;

        // convert plain graph to secure graph.
        edge_block_list.resize(edge_list_size);

        for(size_t i=0; i<edge_list_size; i++){
            aby3::i64Matrix nodes_list(2*l, 1);

            for(size_t j=0; j<l; j++){
                nodes_list(j, 0) = plain_graph.edge_block_list[i][j][0];
                nodes_list(j+l, 0) = plain_graph.edge_block_list[i][j][1];
            }

            // the sbMatrix must be initialized with the correct size.
            edge_block_list[i].resize(2*l, BITSIZE);

            // TODO - optimize
            if(party_info.pIdx == 0){
                party_info.enc->localBinMatrix(*(party_info.runtime), nodes_list, edge_block_list[i]).get();
            }
            else{
                party_info.enc->remoteBinMatrix(*(party_info.runtime), edge_block_list[i]).get();
            }
        }
        return;
    }

    Graph2d(const std::string& meta_data_file, const std::string& edge_block_file, const std::string& node_chunk_file, aby3Info &party_info) : Graph2d(meta_data_file, edge_block_file, party_info) {
        plainGraph2d plain_graph(meta_data_file, edge_block_file, node_chunk_file);

        // encrypt the node list.
        node_chuck_list.resize(b);
        for(size_t i=0; i<b; i++){
            aby3::i64Matrix nodes(k, 1);
            for(size_t j=0; j<k; j++){
                nodes(j, 0) = plain_graph.node_chuck_list[i][j];
            }

            if(party_info.pIdx == 0){
                party_info.enc->localBinMatrix(*(party_info.runtime), nodes, node_chuck_list[i]).get();
            }
            else{
                party_info.enc->remoteBinMatrix(*(party_info.runtime), node_chuck_list[i]).get();
            }
        }

        return;
    }

    void get_node_edges_list(){
        node_edges_list.resize(b);
        for(size_t i=0; i<b; i++){
            node_edges_list[i].resize(b * 2 * l, BITSIZE);
            for(size_t j=0; j<b; j++){
                for(size_t k=0; k<l; k++){
                    // starting node list.
                    node_edges_list[i].mShares[0](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[0](k, 0);
                    node_edges_list[i].mShares[1](j*(l)+k, 0) = edge_block_list[i*b+j].mShares[1](k, 0);
                    // ending node list.
                    node_edges_list[i].mShares[0](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[0](k+l, 0);
                    node_edges_list[i].mShares[1](j*(l)+k+ (b*l), 0) = edge_block_list[i*b+j].mShares[1](k+l, 0);
                }
            }   
        }
        return;
    }

    void check_graph(const std::string& meta_data_file, const std::string& edge_block_file, aby3Info &party_info){
        plainGraph2d plain_graph(meta_data_file, edge_block_file);

        // reveal the secure edge list back to plaintext.
        aby3::i64Matrix start_nodes(edge_list_size * l, 1);
        aby3::i64Matrix end_nodes(edge_list_size * l, 1);

        aby3::sbMatrix start_nodes_sec(edge_list_size * l, 1);
        aby3::sbMatrix end_nodes_sec(edge_list_size * l, 1);

        for(size_t i=0; i<edge_list_size; i++){
            for(size_t j=0; j<l; j++){
                start_nodes_sec.mShares[0](i*l+j, 0) = edge_block_list[i].mShares[0](j, 0);
                start_nodes_sec.mShares[1](i*l+j, 0) = edge_block_list[i].mShares[1](j, 0);
                end_nodes_sec.mShares[0](i*l+j, 0) = edge_block_list[i].mShares[0](j + l, 0);
                end_nodes_sec.mShares[1](i*l+j, 0) = edge_block_list[i].mShares[1](j + l, 0);
            }
        }

        party_info.enc->revealAll(*(party_info.runtime), start_nodes_sec, start_nodes).get();
        party_info.enc->revealAll(*(party_info.runtime), end_nodes_sec, end_nodes).get();
        
        // check whether the secure graph is the same as the plaintext graph.
        bool check_flag = true;
        bool break_outter_loop = false;
        for(size_t i=0; i<edge_list_size; i++){
            for(size_t j=0; j<l; j++){
                if(start_nodes(i*l+j, 0) != plain_graph.edge_block_list[i][j][0] || end_nodes(i*l+j, 0) != plain_graph.edge_block_list[i][j][1]){
                    check_flag = false;
                    break_outter_loop = true;
                    break;
                }
            }
            if(break_outter_loop) break;
        }

        if(party_info.pIdx == 0){
            if(check_flag){
                debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
            }
            else{
                debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");
            }
        }
    }
};

struct GraphAdj {
    size_t v;
    size_t adj_size;
    std::vector<aby3::sbMatrix> adj_list;

    GraphAdj(){}; // default constructor

    GraphAdj(const std::string& meta_data_file, const std::string& data_file, aby3Info &party_info){
        plainGraphAdj plain_graph(meta_data_file, data_file);

        // check the adj graph.
        plain_graph.generate_adj_list();
        v = plain_graph.v;
        adj_size = plain_graph.adj_list.size();
        adj_list.resize(adj_size);

        aby3::i64Matrix full_adj_matrix(adj_size, 1);
        for(size_t i=0; i<adj_size; i++){
            full_adj_matrix(i, 0) = plain_graph.adj_list[i];
        }

        aby3::sbMatrix full_adj_matrix_sec(adj_size, BITSIZE);

        // encrypt the adj graph
        if(party_info.pIdx == 0){
            party_info.enc->localBinMatrix(*(party_info.runtime), full_adj_matrix, full_adj_matrix_sec).get();
        }
        else{
            party_info.enc->remoteBinMatrix(*(party_info.runtime), full_adj_matrix_sec).get();
        }

        for(size_t i=0; i<adj_size; i++){
            adj_list[i].resize(1, BITSIZE);
            adj_list[i].mShares[0](0, 0) = full_adj_matrix_sec.mShares[0](i, 0);
            adj_list[i].mShares[1](0, 0) = full_adj_matrix_sec.mShares[1](i, 0);
        }
        return;
    }

    void check_graph(const std::string& meta_data_file, const std::string& data_file, aby3Info &party_info){
        plainGraphAdj plain_graph(meta_data_file, data_file);
        plain_graph.generate_adj_list();

        // reveal the secure edge list back to plaintext.
        aby3::i64Matrix adj_matrix(plain_graph.adj_list.size(), 1);
        aby3::sbMatrix adj_matrix_sec(plain_graph.adj_list.size(), 1);

        for(size_t i=0; i<plain_graph.adj_list.size(); i++){
            adj_matrix_sec.mShares[0](i, 0) = adj_list[i].mShares[0](0, 0);
            adj_matrix_sec.mShares[1](i, 0) = adj_list[i].mShares[1](0, 0);
        }

        party_info.enc->revealAll(*(party_info.runtime), adj_matrix_sec, adj_matrix).get();
        
        // check whether the secure graph is the same as the plaintext graph.
        bool check_flag = true;
        for(size_t i=0; i<plain_graph.adj_list.size(); i++){
            if(adj_matrix(i, 0) != plain_graph.adj_list[i]){
                check_flag = false;
                break;
            }
        }

        if(party_info.pIdx == 0){
            if(check_flag){
                debug_info("\033[32m The secure graph is the same as the plaintext graph. \033[0m\n");
            }
            else{
                debug_info("\033[31m Error: the secure graph is not the same as the plaintext graph. \033[0m\n");
            }
        }
    }

};

class GraphQueryEngine{
    public:
        Graph2d *graph;
        ABY3SqrtOram *edge_block_oram;
        ABY3SqrtOram *node_edges_oram;
        aby3Info *party_info;

        size_t logb, logk;

        GraphQueryEngine(){}

        GraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& edge_block_file){
            graph = new Graph2d(meta_data_file, edge_block_file, party_info);
            this->party_info = &party_info;

            if(!checkPowerOfTwo(graph->b) || !checkPowerOfTwo(graph->k)){
                THROW_RUNTIME_ERROR("The block size and the chunk size must be power of 2.");
            }

            logb = log2(graph->b);
            logk = log2(graph->k);

            return;
        }

        void edge_block_oram_initialization(const int stash_size, const int pack_size){
            edge_block_oram = new ABY3SqrtOram(graph->edge_list_size, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
            edge_block_oram->initiate(graph->edge_block_list);
            return;
        }

        void node_edges_oram_initialization(const int stash_size, const int pack_size){
            // iniitialize the node edges oram.
            node_edges_oram = new ABY3SqrtOram(graph->b, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

            // construct the node edges data.
            graph->get_node_edges_list();
            node_edges_oram->initiate(graph->node_edges_list);
            return;
        }

        GraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& edge_block_file, const size_t edge_oram_stash_size, const size_t edge_oram_pack_size, const size_t node_oram_stash_size, const size_t node_oram_pack_size) : GraphQueryEngine(party_info, meta_data_file, edge_block_file)
        {
            edge_block_oram_initialization(edge_oram_stash_size, edge_oram_pack_size);
            node_edges_oram_initialization(node_oram_stash_size, node_oram_pack_size);
            return;
        }

        ~GraphQueryEngine(){
            delete graph;
            delete edge_block_oram;
            delete node_edges_oram;
        }

        boolIndex get_block_index(boolIndex node_index){
            boolIndex block_index, left_size;
            bool_shift(party_info->pIdx, node_index, logk, block_index, true); // right shift
            return block_index;
        }

        boolIndex get_edge_block_index(boolIndex starting_node, boolIndex ending_node){
            // this function should be called in plaintext phase.
            // this function is only used for test.
            boolIndex starting_block_index = get_block_index(starting_node);
            boolIndex ending_block_index = get_block_index(ending_node);
            
            boolIndex block_index;
            bool_shift(party_info->pIdx, block_index, logb, starting_block_index, false); // left shift
            aby3::sbMatrix edge_block_index_mat;
            aby3::sbMatrix starting_block_index_mat = starting_block_index.to_matrix();
            aby3::sbMatrix ending_block_index_mat = ending_block_index.to_matrix();
            bool_cipher_add(party_info->pIdx, starting_block_index_mat, ending_block_index_mat, edge_block_index_mat, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

            block_index.from_matrix(edge_block_index_mat);

            return block_index;
        }

        int get_block_index(int node_index){
            return node_index >> logk;
        }

        int get_edge_block_index(int starting_node, int ending_node){
            return ((starting_node >> logk) * graph->b) + (ending_node >> logk);
        }

        aby3::sbMatrix get_edge_block(boolIndex edge_block_idx){
            return edge_block_oram->access(edge_block_idx);
        }

        aby3::sbMatrix get_node_edges(boolIndex node_idx){
            return node_edges_oram->access(node_idx);
        }

        void print_configs(std::ostream& stream){
            stream << "v : " << graph->v << std::endl;
            stream << "e : " << graph->e << std::endl;
            stream << "b : " << graph->b << std::endl;
            stream << "l : " << graph->l << std::endl;
            stream << "k : " << graph->k << std::endl;
            stream << "se : " << edge_block_oram->S << std::endl;
            stream << "pe : " << edge_block_oram->pack << std::endl;
            stream << "sn : " << node_edges_oram->S << std::endl;
            stream << "pn : " << node_edges_oram->pack << std::endl;
        }
};

class AdjGraphQueryEngine{

public:
    GraphAdj *graph;
    ABY3SqrtOram *edge_oram;
    ABY3SqrtOram *node_oram;
    aby3Info *party_info;

    size_t v;

    AdjGraphQueryEngine(){}

    AdjGraphQueryEngine(aby3Info &party_info, const std::string& meta_data_file, const std::string& data_file){
        graph = new GraphAdj(meta_data_file, data_file, party_info);
        this->party_info = &party_info;
        v = graph->v;

        if(!checkPowerOfTwo(v)){ // as we only support 2^m indexed ORAM for efficiency.
            THROW_RUNTIME_ERROR("The number of nodes must be power of 2.");
        }
        return;
    }

    void edge_oram_initialization(const int stash_size, const int pack_size){
        edge_oram = new ABY3SqrtOram(graph->adj_size, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
        edge_oram->initiate(graph->adj_list);
        return;
    }

    void node_oram_initialization(const int stash_size, const int pack_size){

        // organize the node-edges data structure.
        std::vector<aby3::sbMatrix> node_edges_list(v);

        for(size_t i=0; i<v; i++){
            aby3::sbMatrix node_edges(v, BITSIZE);
            for(size_t j=0; j<v; j++){
                node_edges.mShares[0](j, 0) = graph->adj_list[i*v+j].mShares[0](0, 0);
                node_edges.mShares[1](j, 0) = graph->adj_list[i*v+j].mShares[1](0, 0);
            }
            node_edges_list[i] = node_edges;
        }

        node_oram = new ABY3SqrtOram(v, stash_size, pack_size, party_info->pIdx, *(party_info->enc), *(party_info->eval), *(party_info->runtime));
        node_oram->initiate(node_edges_list);

        return;
    }

    aby3::sbMatrix get_target_edge(boolIndex edge_index){
        return edge_oram->access(edge_index);
    }

    aby3::sbMatrix get_target_node(boolIndex node_index){
        return node_oram->access(node_index);
    }

    boolIndex get_logical_edge_index(boolIndex start_node, boolIndex end_node){
        boolIndex logical_edge_index;
        bool_shift(party_info->pIdx, start_node, log2(v), logical_edge_index, false); // left shift
        aby3::sbMatrix edge_index_mat;
        aby3::sbMatrix starting_index_mat = logical_edge_index.to_matrix();
        aby3::sbMatrix ending_index_mat = end_node.to_matrix();
        bool_cipher_add(party_info->pIdx, starting_index_mat, ending_index_mat, edge_index_mat, *(party_info->enc), *(party_info->eval), *(party_info->runtime));

        logical_edge_index.from_matrix(edge_index_mat);

        return logical_edge_index;
    
    }

};

aby3::sbMatrix get_target_node_mask(boolIndex target_start_node, aby3::sbMatrix& node_block, aby3Info &party_info);

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,
                         boolIndex logical_edge_block_index,
                         GraphQueryEngine &GQEngine);

aby3::si64Matrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,AdjGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_edge_count(boolIndex boolIndex, AdjGraphQueryEngine &GQEngine);

aby3::sbMatrix outting_neighbors(boolIndex node_index,AdjGraphQueryEngine &GQEngine);