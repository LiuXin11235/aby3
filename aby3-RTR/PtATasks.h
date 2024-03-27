#include "GeneralPTA.h"
#include "DataGenerator.h"

class CipherIndex : SubABY3Task<aby3::si64, int, aby3::si64, aby3::si64> {
    using SubABY3Task<aby3::si64, int, aby3::si64, aby3::si64>::SubABY3Task;

    CipherIndex(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval)
        : SubABY3Task<aby3::si64, int, aby3::si64, aby3::si64>(optimal_block, task_id, pIdx, enc, runtime, eval) {
        this->have_selective = true;
    }

    virtual void compute_local_table(std::vector<aby3::si64>& expandX, std::vector<int>& expandY, std::vector<aby3::si64>& local_table, BlockInfo* binfo) override {
        aby3::u64 block_length = binfo->block_len;
        
        // pairwise comparison
        aby3::sbMatrix partTable;
        vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *this->eval, *this->runtime, *this->enc);

        // pairwise abmul.
        aby3::si64Matrix expandV(block_length, 1);
        for (size_t i = 0; i < block_length; i++) expandV(i, 0, this->selectV[binfo->t_start + i]);
        cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *this->eval, *this->enc, *this->runtime);

        //trans to local_table.
        for (size_t i = 0; i < block_length; i++) local_table[i] = expandV(i, 0);
    }

    virtual void partical_reduction(std::vector<aby3::si64>& resLeft,
                                  std::vector<aby3::si64>& resRight,
                                  std::vector<aby3::si64>& local_res,
                                  BlockInfo* binfo) override {
        for (size_t i = 0; i < resLeft.size(); i++) local_res[i] = resLeft[i] + resRight[i];
    }

    std::tuple<std::vector<aby3::si64>, std::vector<int>, std::vector<aby3::si64>> data_loading(std::string data_folder){
        // load the input data.
        raise NotImplementedError("data_loading through file");
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