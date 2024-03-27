#include "GeneralPTA.h"
#include "DataGenerator.h"

const aby3::si64 GET_ZERO_SHARE = [] {
    aby3::si64 dval;
    dval.mData[0] = 0;
    dval.mData[1] = 0;
    return dval;
}();

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
        aby3::si64Matrix expandV(block_length, 1);
        for (size_t i = 0; i < block_length; i++) expandV(i, 0, this->selectV[binfo->t_start + i]);
        cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *this->eval, *this->enc, *this->runtime);

        //trans to local_table.
        for (size_t i = 0; i < block_length; i++) local_table[i] = expandV(i, 0);
    }

    virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
        for (size_t i = 0; i < resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
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
