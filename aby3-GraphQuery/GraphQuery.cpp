#include "Graph.h"

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

aby3::sbMatrix outting_edge_count(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){

    // get the node block from the GraohQueryEngine, which size is b * 2l.
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

    // aggregate the eq_res for the final result.
    eq_res.resize(eq_res.rows(), BITSIZE);
    aby3::sbMatrix res(1, BITSIZE);

    bool_aggregation(GQEngine.party_info->pIdx, eq_res, res, *(GQEngine.party_info->enc), *(GQEngine.party_info->eval), *(GQEngine.party_info->runtime), "ADD");

    return res;
}

// TODO: neighbors finding query.
aby3::sbMatrix outting_neighbors(boolIndex node_index, boolIndex logical_node_block_index, GraphQueryEngine &GQEngine){
    THROW_RUNTIME_ERROR("Not implemented yet.");

    return aby3::sbMatrix(1, 1);
}