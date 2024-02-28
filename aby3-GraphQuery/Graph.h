#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include "pGraph.h"
#include "../aby3-Basic/SqrtOram.h"

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
};

boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,
                         boolIndex logical_edge_block_index,
                         GraphQueryEngine &GQEngine);

aby3::sbMatrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine);