#pragma once
#include "GeneralPTA.h"
#include "DataGenerator.h"

const aby3::si64 GET_ZERO_SHARE = [] {
    aby3::si64 dval;
    dval.mData[0] = 0;
    dval.mData[1] = 0;
    return dval;
}();

template<typename T>
struct indexData{
  T value;
  aby3::si64 index;
};

// functional class override.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class CipherIndex : public SubTask<NUMX, NUMY, NUMT, NUMR> { 

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    CipherIndex(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = true;
    }

    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        // pairwise comparison
        aby3::sbMatrix partTable(block_length, 1); 
        vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *(this->eval), *(this->runtime));

        // pairwise abmul.
        aby3::si64Matrix expandV(block_length, 1), p_v_after_mul(block_length, 1);
        for(size_t i=0; i<block_length; i++) {
            expandV.mShares[0](i, 0) = this->selectV[binfo->t_start + i].mData[0];
            expandV.mShares[1](i, 0) = this->selectV[binfo->t_start + i].mData[1];
        }

        cipher_mul(this->pIdx, expandV, partTable, p_v_after_mul, *this->eval, *this->enc, *this->runtime);

        //trans to local_table.
        for (size_t i = 0; i < block_length; i++) {
            local_table[i].mData[0] = p_v_after_mul.mShares[0](i, 0);
            local_table[i].mData[1] = p_v_after_mul.mShares[1](i, 0);
        }
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for (size_t i = 0; i < resLeft.size(); i++) {
            local_res[i] = resLeft[i] + resRight[i];
        }
    }

    std::tuple<std::vector<aby3::si64>, std::vector<int>, std::vector<aby3::si64>> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::tuple<std::vector<aby3::si64>, std::vector<int>, std::vector<aby3::si64>> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecX = generate_vector_si64(this->n, this->pIdx, *this->enc, *this->runtime);
        std::vector<int> vecY = generate_vector_int(partial_len);
        std::vector<aby3::si64> vecV = generate_vector_si64(partial_len,  this->pIdx, *this->enc, *this->runtime);
        return std::make_tuple(vecX, vecY, vecV);
    }
};

// Max function class
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Max : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Max(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        local_table = expandY;
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        aby3::sbMatrix comp_res;
        vector_cipher_gt(this->pIdx, resLeft, resRight, comp_res, *(this->eval), *(this->enc), *(this->runtime));
        aby3::si64Matrix mul_mat(resLeft.size(), 1), mul_res(resLeft.size(), 1);
        for(size_t i=0; i<resLeft.size(); i++) mul_mat(i, 0, resRight[i] - resLeft[i]);
        cipher_mul_seq(this->pIdx, mul_mat, comp_res, mul_res, *(this->eval), *(this->enc), *(this->runtime));
        for(size_t i=0; i<resLeft.size(); i++) local_res[i] = resLeft[i] + mul_res(i, 0);
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::vector<aby3::si64> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return vecY;
    }

};

// Rank function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Rank : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Rank(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        aby3::sbMatrix partTable(block_length, 1);
        vector_cipher_gt(this->pIdx, expandX, expandY, local_table, *(this->eval), *(this->enc), *(this->runtime));
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for(size_t i=0; i<resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::pair<std::vector<aby3::si64>, std::vector<aby3::si64>> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecX = generate_vector_si64(this->n, this->pIdx, *this->enc, *this->runtime);
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return std::make_pair(vecX, vecY);
    }

};

// Sum function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class Sum : public SubTask<NUMX, NUMY, NUMT, NUMR> {

public:
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    Sum(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        local_table = expandY;
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for(size_t i=0; i<resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
        return;
    }

    std::vector<aby3::si64> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::vector<aby3::si64> data_loading(){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<aby3::si64> vecY = generate_vector_si64(partial_len, this->pIdx, *this->enc, *this->runtime);
        return vecY;
    }
};

// BioMetric function class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class BioMetric : public SubTask<NUMX, NUMY, NUMT, NUMR> {
    //aby3 info
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    BioMetric(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval) : 
        pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
            this->have_selective = false;
    }

    // functional code.
    virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        size_t k = expandX[0].size();

        // flat the two-dimensional inputs.
        size_t exp_len = expandX.size()*k;
        std::vector<typename NUMX::value_type> flatX, flatY;
        for (const auto& innerVec : expandX){
            for (const auto& element : innerVec) flatX.push_back(element);
        }
        for (const auto& innerVec : expandY){
            for (const auto& element : innerVec) flatY.push_back(element);
        }
        // vector mul
        if(std::is_same<typename NUMX::value_type, aby3::si64>::value){
            vector_mean_square(this->pIdx, flatX, flatY, flatX, *(this->eval), *(this->enc), *(this->runtime));
        }

        // enc index.
        aby3::i64Matrix indexMat(expandX.size(), 1);
        aby3::si64Matrix sindex(expandX.size(), 1);
        for(int i=0; i<expandX.size(); i++) indexMat(i, 0) = binfo->t_start + i;
        if(this->pIdx == 0){
            this->enc->localIntMatrix(*(this->runtime), indexMat, sindex).get();
        }
        else{
            this->enc->remoteIntMatrix(*(this->runtime), sindex).get();
        }
        
        // reduce to local table
        for(int i=0; i<expandX.size(); i++){
            local_table[i].value = this->initial_value.value;
            local_table[i].index = sindex(i, 0);
            for (int j=0; j<k; j++){
                local_table[i].value = local_table[i].value + flatX[j];
            }
        }
        return;
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        aby3::sbMatrix comp_res;
        std::vector<aby3::si64> resLeft_value(resLeft.size());
        std::vector<aby3::si64> resRight_value(resRight.size());
        // initiate the value
        for(int i=0; i<resLeft.size(); i++){
        resLeft_value[i] = resLeft[i].value; resRight_value[i] = resRight[i].value;
        }
        vector_cipher_gt(this->pIdx, resLeft_value, resRight_value, comp_res, *(this->eval), *(this->enc), *(this->runtime));

        // multiply for value extraction.
        aby3::si64Matrix matSub;
        matSub.resize(resLeft.size(), 1);
        for(int i=0; i<resLeft.size(); i++) matSub(static_cast<aby3::u64>(i), static_cast<aby3::u64>(0), resRight[i].value - resLeft[i].value);
        cipher_mul_seq(this->pIdx, matSub, comp_res, matSub, *(this->eval), *(this->enc), *(this->runtime));

        aby3::si64Matrix matIndex;
        matIndex.resize(resLeft.size(), 1);
        for(int i=0; i<matIndex.size(); i++) matIndex(i, 0, resRight[i].index - resLeft[i].index);
        cipher_mul_seq(this->pIdx, matIndex, comp_res, matIndex, *(this->eval), *(this->enc), *(this->runtime));

        // compute the final result
        for(int i=0; i<resLeft.size(); i++){
            local_res[i].value = resRight[i].value - matSub(i, 0);
            local_res[i].index = resRight[i].index - matIndex(i, 0);
        }
        return;
    }

    std::pair<std::vector<std::vector<aby3::si64>>, std::vector<std::vector<aby3::si64>>> data_loading(std::string data_folder){
        // load the input data.
        THROW_RUNTIME_ERROR("data_loading through file is not implemented yet.");
    }

    std::pair<std::vector<std::vector<aby3::si64>>, std::vector<std::vector<aby3::si64>>> data_loading(size_t k=1){
        // directly generate the input data.
        size_t partial_len = this->get_partial_m_lens();
        std::vector<std::vector<aby3::si64>> vecX = generate_vector_si64(this->n * this->k, this->pIdx, *this->enc, *this->runtime);
        std::vector<std::vector<aby3::si64>> vecY = generate_vector_si64(partial_len * this->k, this->pIdx, *this->enc, *this->runtime);

        std::vector<std::vector<aby3::si64>> vecX2d(this->n, std::vector<aby3::si64>(k));
        std::vector<std::vector<aby3::si64>> vecY2d(partial_len, std::vector<aby3::si64>(k));
        for(size_t i=0; i<this->n; i++){
            for(size_t j=0; j<k; j++){
                vecX2d[i][j] = vecX[i*k + j];
            }
        }
        for(size_t i=0; i<partial_len; i++){
            for(size_t j=0; j<k; j++){
                vecY2d[i][j] = vecY[i*k + j];
            }
        }
        return std::make_pair(vecX2d, vecY2d);
    }
};