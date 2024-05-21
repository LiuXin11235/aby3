#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <cassert>
#include <string>

struct plainGraph2d{
    size_t v, e;
    size_t b, k, l;
    size_t edge_list_size;

    std::vector<std::vector<int>> node_chuck_list;
    std::vector<std::vector<std::array<int, 2>>> edge_block_list;

    plainGraph2d(){}; // default constructor

    plainGraph2d(const std::string& meta_data_file, const std::string& edge_block_file){
        // load the graph meta data.
        std::ifstream meta(meta_data_file);
        meta >> v >> e >> b >> k >> l;
        this->edge_list_size = b * b;
        node_chuck_list.resize(b);
        edge_block_list.resize(edge_list_size);

        // load the edges.
        std::ifstream edge_block(edge_block_file);
        for (size_t i = 0; i < edge_list_size; i++) {
            edge_block_list[i].resize(l);
            for (size_t j = 0; j < l; j++) {
                edge_block >> edge_block_list[i][j][0] >> edge_block_list[i][j][1];
            }
        }
        return;
    }
    
    plainGraph2d(const std::string& meta_data_file, const std::string& edge_block_file, const std::string& node_chunk_file) : plainGraph2d(meta_data_file, edge_block_file)
    {
        // load the nodes.
        std::ifstream node_chunk(node_chunk_file);
        for (size_t i = 0; i < b; i++) {
            node_chuck_list[i].resize(k);
            for (size_t j = 0; j < k; j++) {
                node_chunk >> node_chuck_list[i][j];
            }
        }
        return;
    }

    void printGraphMeta(){
        std::cout << "v: " << this->v << std::endl;
        std::cout << "e: " << this->e << std::endl;
        std::cout << "b: " << this->b << std::endl;
        std::cout << "k: " << this->k << std::endl;
        std::cout << "l: " << this->l << std::endl;
        std::cout << "edge_list_size: " << this->edge_list_size << std::endl;
    }

    void printEdgeList(){
        for(size_t i=0; i<edge_list_size; i++){
            // std::cout << "edge_block_list[" << i << "]: " << std::endl;
            for(size_t j=0; j<l; j++){
                std::cout << edge_block_list[i][j][0] << " " << edge_block_list[i][j][1] << std::endl;
            }
        }
    }

    std::vector<int> get_edge_block(int starting_node, int ending_node, bool block_flag=true){
        int starting_chunk, ending_chunk;
        if(!block_flag){
            starting_chunk = starting_node / k;
            ending_chunk = ending_node / k;
        }
        else{
            starting_chunk = starting_node;
            ending_chunk = ending_node;
        }

        std::vector<int> edge_block(2*l);
        for(size_t i=0; i<l; i++){
            edge_block[i] = edge_block_list[starting_chunk * b + ending_chunk][i][0];
            edge_block[i + l] = edge_block_list[starting_chunk * b + ending_chunk][i][1];
        }

        return edge_block;
    }

    std::vector<int> get_node_chunk(int starting_node, bool block_flag=true){
        int starting_chunk;
        if(!block_flag){
            starting_chunk = starting_node / k;
        }
        else{
            starting_chunk = starting_node;
        }

        std::vector<int> node_chunk(b * 2 * l);
        for(size_t i=0; i<b; i++){
            std::vector<int> edge_block = get_edge_block(starting_chunk, i);
            for(size_t j=0; j<l; j++){
                node_chunk[i * (l) + j] = edge_block[j]; // starting nodes.
                node_chunk[i * (l) + j + (b*l)] = edge_block[j + l]; // ending nodes.
            }
        }

        return node_chunk;
    }

    std::vector<int> get_outting_neighbors(int starting_node){
        std::vector<int> neighbors;
        std::vector<int> node_chunk = get_node_chunk(starting_node);
        for(size_t i=0; i<b*l; i++){
            if(node_chunk[i] == starting_node){
                neighbors.push_back(node_chunk[b*l + i]);
            }
        }
        return neighbors;
    }
};

struct plainGraphAdj{
    size_t v;
    size_t unique_edges = 0;
    std::vector<int> starting_node_list;
    std::vector<int> ending_node_list;  
    std::vector<int> adj_list;

    plainGraphAdj(){}; // default constructor

    plainGraphAdj(const std::string& meta_data_file, const std::string& edge_list_file){
        // load the graph meta data.
        std::ifstream meta(meta_data_file);
        meta >> v;
        starting_node_list.resize(v);
        ending_node_list.resize(v);

        // load the edges.
        std::ifstream edge_list(edge_list_file);
        for (size_t i = 0; i < v; i++) {
            edge_list >> starting_node_list[i] >> ending_node_list[i];
        }
        return;
    }

    void generate_adj_list(){
        assert(!starting_node_list.empty() && !ending_node_list.empty());
        assert(starting_node_list.size() == ending_node_list.size());

        adj_list.resize(v*v, 0);
        for(size_t i=0; i<starting_node_list.size(); i++){
            int s = starting_node_list[i], e = ending_node_list[i];
            if(adj_list[s*v + e] == 0){
                unique_edges++;
            }
            adj_list[s*v + e] += 1;
        }
        return;
    }

    bool edge_existence(int starting_node, int ending_node){
        return adj_list[starting_node*v + ending_node] > 0;
    }

    int outting_neighbors_count(int starting_node){
        int count = 0;
        for(size_t i=0; i<v; i++){
            count += adj_list[starting_node*v + i];
        }
        return count;
    }

    std::vector<int> get_outting_neighbors(int starting_node){
        std::vector<int> neighbors(v, 0);
        for(size_t i=0; i<v; i++){
            if(adj_list[starting_node*v + i] > 0){
                neighbors[i] = 1;
            }
        }
        return neighbors;
    }

    void printGraphMeta(){
        std::cout << "v: " << this->v << std::endl;
        std::cout << "unique_edges: " << this->unique_edges << std::endl;
    }
};