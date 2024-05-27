#include "Graph.h"

aby3::sbMatrix get_target_node_mask(boolIndex target_start_node, aby3::sbMatrix& node_block, aby3Info &party_info){

    // check whether node_block contains even rows.
    if(checkEven(node_block.rows()) == false){
        THROW_RUNTIME_ERROR("The node block should contain even rows.");
    }

    size_t edge_num = node_block.rows() / 2;

    // process the node block for the final result.
    aby3::sbMatrix expanded_node_mat(edge_num, BITSIZE);
    aby3::sbMatrix starting_node_mat(edge_num, BITSIZE);

    for(int i=0; i<2; i++){ // iterate for the two shares.
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + edge_num, expanded_node_mat.mShares[i].begin());
        for(size_t j=0; j<edge_num; j++){
            starting_node_mat.mShares[i](j, 0) = target_start_node.indexShares[i];
        }
    }

    // comparing the starting nodes with the target starting node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(party_info.pIdx, expanded_node_mat, starting_node_mat, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    return eq_res;
}

aby3::sbMatrix get_unique_ending_nodes_in_edge_list(boolIndex target_start_node, aby3::sbMatrix& starting_nodes, aby3::sbMatrix& ending_nodes, aby3Info& party_info){

    size_t edge_num = starting_nodes.rows();    
    aby3::sbMatrix target_starting_nodes(edge_num, BITSIZE);
    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            target_starting_nodes.mShares[j](i, 0) = target_start_node.indexShares[j];
        }
    }

    // get the mask indicating which elements are the target node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(party_info.pIdx, starting_nodes, target_starting_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // expand the length to match with the node id.
    aby3::sbMatrix aligned_target_node_mat(edge_num, BITSIZE);
    for(int i=0; i<2; i++){
        for(size_t j=0; j<edge_num; j++){
            aligned_target_node_mat.mShares[i](j, 0) = (eq_res.mShares[i](j, 0) == 1) ? -1 : 0;
        }
    }

    // multiply with the real starting_nodes to extract the target nodes.
    bool_cipher_and(party_info.pIdx, aligned_target_node_mat, ending_nodes, aligned_target_node_mat, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // sort the nodes for grouping.
    std::vector<aby3::sbMatrix> masked_nodes(edge_num);
    for(size_t i=0; i<edge_num; i++){
        masked_nodes[i].resize(1, BITSIZE);
        for(int j=0; j<2; j++){
            masked_nodes[i].mShares[j](0, 0) = aligned_target_node_mat.mShares[j](i, 0);
        }
    }

    quick_sort(masked_nodes, party_info.pIdx, *(party_info.enc), *(party_info.eval), *(party_info.runtime), 32);

    // differente EQ to filter out the fliping nodes.
    aby3::sbMatrix left_shifted_nodes(edge_num-1, BITSIZE);
    aby3::sbMatrix right_shifted_nodes(edge_num-1, BITSIZE);

    for (size_t i = 0; i < edge_num-1; i++)
    {
        left_shifted_nodes.mShares[0](i, 0) = masked_nodes[i].mShares[0](0, 0);
        left_shifted_nodes.mShares[1](i, 0) = masked_nodes[i].mShares[1](0, 0);
        right_shifted_nodes.mShares[0](i, 0) = masked_nodes[i+1].mShares[0](0, 0);
        right_shifted_nodes.mShares[1](i, 0) = masked_nodes[i+1].mShares[1](0, 0);
    }

    bool_cipher_eq(party_info.pIdx, left_shifted_nodes, right_shifted_nodes, eq_res, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    eq_res.resize(eq_res.rows(), BITSIZE);
    for(size_t i=0; i<eq_res.rows(); i++){
        for(int j=0; j<2; j++){
            eq_res.mShares[j](i, 0) = (eq_res.mShares[j](i, 0) == 1) ? -1 : 0;
        }
    }
    eq_res.resize(eq_res.rows()+1, BITSIZE);
    boolIndex true_share(0, party_info.pIdx);
    eq_res.mShares[0](eq_res.rows()-1, 0) = true_share.indexShares[0];
    eq_res.mShares[1](eq_res.rows()-1, 0) = true_share.indexShares[1];

    bool_cipher_not(party_info.pIdx, eq_res, eq_res);

    for(size_t i=0; i<edge_num; i++){
        for(int j=0; j<2; j++){
            aligned_target_node_mat.mShares[j](i, 0) = masked_nodes[i].mShares[j](0, 0);
        }
    }

    aby3::sbMatrix filtered_nodes(edge_num, BITSIZE);
    bool_cipher_and(party_info.pIdx, eq_res, aligned_target_node_mat, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    // shuffle the masks and target nodes for privacy.
    efficient_shuffle(filtered_nodes, party_info.pIdx, filtered_nodes, *(party_info.enc), *(party_info.eval), *(party_info.runtime));

    return filtered_nodes;
}

// functions using GraphQueryEngine (based on Graph2D).
boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,
                         boolIndex logical_edge_block_index,
                         GraphQueryEngine &GQEngine) {
    // get the edge block from the GraphQueryEngine, which size is 2*l,
    // containing l beginning node and l ending node.
    aby3::sbMatrix edge_block =
        GQEngine.get_edge_block(logical_edge_block_index);

    // process the edge block for the final result.
    aby3::sbMatrix expanded_node_mat(GQEngine.graph->l * 2, BITSIZE);
    for (size_t i = 0; i < GQEngine.graph->l; i++) {
        expanded_node_mat.mShares[0](i, 0) = starting_node.indexShares[0];
        expanded_node_mat.mShares[1](i, 0) = starting_node.indexShares[1];
        expanded_node_mat.mShares[0](i + GQEngine.graph->l, 0) =
            ending_node.indexShares[0];
        expanded_node_mat.mShares[1](i + GQEngine.graph->l, 0) =
            ending_node.indexShares[1];
    }

    // comparing the starting nodes and ending nodes with the target starting
    // node and rge ending node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, edge_block, expanded_node_mat,
                   eq_res, *(GQEngine.party_info->enc),
                   *(GQEngine.party_info->eval),
                   *(GQEngine.party_info->runtime));

    // check whether there exist pair (starting_node, ending_node) in the edge
    // block.
    aby3::sbMatrix starting_match(GQEngine.graph->l, 1),
        ending_match(GQEngine.graph->l, 1);

    for (int i = 0; i < 2; i++) {
        std::copy(eq_res.mShares[i].begin(),
                  eq_res.mShares[i].begin() + GQEngine.graph->l,
                  starting_match.mShares[i].begin());
        std::copy(eq_res.mShares[i].begin() + GQEngine.graph->l,
                  eq_res.mShares[i].end(), ending_match.mShares[i].begin());
    }

    aby3::sbMatrix matching_res(GQEngine.graph->l, 1);
    bool_cipher_and(GQEngine.party_info->pIdx, starting_match, ending_match,
                    matching_res, *(GQEngine.party_info->enc),
                    *(GQEngine.party_info->eval),
                    *(GQEngine.party_info->runtime));

    // aggregate for the final result, all using log-round OP computations.
    aby3::sbMatrix res(1, 1);
    bool_aggregation(GQEngine.party_info->pIdx, matching_res, res,
                        *(GQEngine.party_info->enc),
                        *(GQEngine.party_info->eval),
                        *(GQEngine.party_info->runtime), "OR");

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    return query_res;
}

aby3::si64Matrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){

    // get the node block from the GraphQueryEngine, which size is b * 2l.
    aby3::sbMatrix node_block = GQEngine.get_node_edges(logical_node_block_index);

    // process the node block for the final result.
    aby3::sbMatrix expanded_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);
    aby3::sbMatrix starting_node_mat(GQEngine.graph->l * GQEngine.graph->b, BITSIZE);

    // TODO: node chunk data structure.
    for(int i=0; i<2; i++){
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + GQEngine.graph->l * GQEngine.graph->b, expanded_node_mat.mShares[i].begin());
        for(size_t j=0; j<GQEngine.graph->l * GQEngine.graph->b; j++){
            starting_node_mat.mShares[i](j, 0) = node_index.indexShares[i];
        }
    }

    // comparing the starting nodes with the target starting node.
    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, expanded_node_mat, starting_node_mat, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    // trans2 A shares.
    aby3::si64Matrix eq_a_res(eq_res.rows(), 1);
    bool2arith(GQEngine.party_info->pIdx, eq_res, eq_a_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::si64Matrix res(1, 1);
    arith_aggregation(GQEngine.party_info->pIdx, eq_a_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}


aby3::sbMatrix outting_neighbors(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){

    // get the node block from the GraphQueryEngine, which size is b * 2l.
    aby3::sbMatrix node_block = GQEngine.get_node_edges(logical_node_block_index);
    size_t edge_num = GQEngine.graph->l * GQEngine.graph->b;

    aby3::sbMatrix starting_nodes(edge_num, BITSIZE);
    aby3::sbMatrix ending_nodes(edge_num, BITSIZE);

    for(int i=0; i<2; i++){
        std::copy(node_block.mShares[i].begin(), node_block.mShares[i].begin() + edge_num, starting_nodes.mShares[i].begin());
        std::copy(node_block.mShares[i].begin() + edge_num, node_block.mShares[i].end(), ending_nodes.mShares[i].begin());
    }

    return get_unique_ending_nodes_in_edge_list(node_index, starting_nodes, ending_nodes, *(GQEngine.party_info));
}

// functions using AdjGraphQueryEngine.
boolShare edge_existance(boolIndex starting_node, boolIndex ending_node,AdjGraphQueryEngine &GQEngine){
    boolIndex logical_edge_index = GQEngine.get_logical_edge_index(starting_node, ending_node);

    aby3::sbMatrix edge_info = GQEngine.get_target_edge(logical_edge_index);

    // whdther edge_info > 0 or not.
    aby3::sbMatrix res(1, 1);
    aby3::sbMatrix zero_share(edge_info.rows(), edge_info.bitCount());
    for(int i=0; i<edge_info.rows(); i++){
        zero_share.mShares[0](i, 0) = 0;
        zero_share.mShares[1](i, 0) = 0;
    }

    bool_cipher_lt(GQEngine.party_info->pIdx, zero_share, edge_info, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    return query_res;
}

aby3::sbMatrix outting_edge_count(boolIndex node_index, AdjGraphQueryEngine &GQEngine){
    aby3::sbMatrix res(1, 1);
    aby3::sbMatrix zero_share(1, 1);
    zero_share.mShares[0](0, 0) = 0;
    zero_share.mShares[1](0, 0) = 0;

    aby3::sbMatrix node_info = GQEngine.get_target_node(node_index);

    bool_aggregation(GQEngine.party_info->pIdx, node_info, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

aby3::sbMatrix outting_neighbors(boolIndex node_index, AdjGraphQueryEngine &GQEngine){
    aby3::sbMatrix res = GQEngine.get_target_node(node_index);

    // filterout the > 0 elements all to 1.
    aby3::sbMatrix zero_share(res.rows(), res.bitCount());
    for(int i=0; i<res.rows(); i++){
        zero_share.mShares[0](i, 0) = 0;
        zero_share.mShares[1](i, 0) = 0;
    }
    bool_cipher_lt(GQEngine.party_info->pIdx, zero_share, res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    return res;
}


boolShare edge_existance(boolIndex starting_node, boolIndex ending_node, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    aby3::sbMatrix expand_ending_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = starting_node.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = starting_node.indexShares[1];
        expand_ending_node.mShares[0](i, 0) = ending_node.indexShares[0];
        expand_ending_node.mShares[1](i, 0) = ending_node.indexShares[1];
    }

    aby3::sbMatrix eq_res_starts, eq_res_ends;
    bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.starting_node_list, expand_starting_node, eq_res_starts, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.ending_node_list, expand_ending_node, eq_res_ends, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::sbMatrix full_comp_res;
    bool_cipher_and(GQEngine.party_info->pIdx, eq_res_starts, eq_res_ends, full_comp_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));

    aby3::sbMatrix res(1, 1);
    bool_aggregation(GQEngine.party_info->pIdx, full_comp_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "OR");

    boolShare query_res;
    query_res.from_matrix(res.mShares[0](0, 0), res.mShares[1](0, 0));

    return query_res;
}

aby3::sbMatrix outting_edge_count(boolIndex boolIndex, ListGraphQueryEngine &GQEngine){
    aby3::sbMatrix expand_starting_node(GQEngine.e, BITSIZE);
    for(size_t i=0; i<GQEngine.e; i++){
        expand_starting_node.mShares[0](i, 0) = boolIndex.indexShares[0];
        expand_starting_node.mShares[1](i, 0) = boolIndex.indexShares[1];
    }

    aby3::sbMatrix eq_res;
    bool_cipher_eq(GQEngine.party_info->pIdx, GQEngine.starting_node_list, expand_starting_node, eq_res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime));
    eq_res.resize(eq_res.rows(), BITSIZE);
    for(size_t i=0; i<eq_res.rows(); i++){
        for(int j=0; j<2; j++){
            eq_res.mShares[j](i, 0) = (eq_res.mShares[j](i, 0) == 1) ? 1 : 0;
        }
    }

    aby3::sbMatrix res(1, BITSIZE);
    bool_aggregation(GQEngine.party_info->pIdx, eq_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

aby3::sbMatrix outting_neighbors(boolIndex node_index, ListGraphQueryEngine &GQEngine){

    return get_unique_ending_nodes_in_edge_list(node_index, GQEngine.starting_node_list, GQEngine.ending_node_list, *(GQEngine.party_info));
}