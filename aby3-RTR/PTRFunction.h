#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Network/IOService.h>
#include <mpi.h>

#include "BuildingBlocks.h"
#include "./Pair_then_Reduce/include/datatype.h"
#include "./Pair_then_Reduce/include/tasks.h"

#define OPTIMAL_BLOCK 78
#define TASKS 500

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class SubIndex : public SubTask<NUMX, NUMY, NUMT, NUMR>{

public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor* enc;
  aby3::Sh3Runtime* runtime;
  aby3::Sh3Evaluator* eval;

  using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

  SubIndex(const size_t optimal_block, const int task_id, const int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime& runtime, aby3::Sh3Evaluator& eval):
    pIdx(pIdx), enc(&enc), runtime(&runtime), eval(&eval),
    SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id){
        this->free_combine = true;
        this->look_ahead = 0;
        this->have_selective = true;
    }

protected:

  virtual void compute_local_table(std::vector<NUMX>& expandX, std::vector<NUMY>& expandY, std::vector<NUMT>& local_table, BlockInfo* binfo){
    aby3::u64 block_length = binfo->block_len;
    // pairwise vector equal test.
    aby3::sbMatrix partTable;
    vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *(this->eval), *(this->runtime));

    // pairwise vector abmul.
    aby3::si64Matrix expandV;
    expandV.resize(block_length, 1);
    for(size_t i=0; i<block_length; i++){
      expandV(i, 0, this->selectV[binfo->t_start+i]);
    }
    cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *(this->eval), *(this->enc), *(this->runtime));

    // transfor to the local_table
    for(size_t i=0; i<block_length; i++){
      local_table[i] = expandV(i, 0);
    }
  }

  virtual void partical_reduction(std::vector<NUMT>& local_table, std::vector<NUMR>& local_res, BlockInfo* binfo) override {
    for(int i=0; i<binfo->res_len; i++){
            size_t row_start = binfo->table_rows[i], row_end = binfo->table_rows[i+1];
            local_res[i] = local_table[row_start];
            for(int j=row_start+1; j<row_end; j++){
                local_res[i] = local_res[i] + local_table[j];
            }
        }
    return;
  }

  virtual void direct_combine(std::vector<NUMR>& local_res, BlockInfo* binfo) override {
    for(int i=0; i<binfo->res_len; i++){
            this->partical_res[i+binfo->r_start] = this->partical_res[i+binfo->r_start] + local_res[i];
        }
        return;
  }

  virtual void final_combine() override {}

};


template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,   
          template<typename, typename, typename, typename> class TASK>
class SecretIndex : public PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>{

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor& enc;
    aby3::Sh3Runtime& runtime;
    aby3::Sh3Evaluator& eval;

    using PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::PTRTask;
    SecretIndex(int tasks, size_t optimal_block, const int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime& runtime, aby3::Sh3Evaluator& eval):
    pIdx(pIdx), enc(enc), runtime(runtime), eval(eval),
    PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block){}

    SecretIndex(const int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime& runtime, aby3::Sh3Evaluator& eval, std::string deployment_profile_name="./Config/profile.json"):
    pIdx(pIdx), enc(enc), runtime(runtime), eval(eval),
    PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(deployment_profile_name){}

    using PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::task_split;  
    void task_split(size_t table_size, size_t task_length) override {
      for(int i=0; i<this->total_tasks; i++){
          size_t table_start = i*task_length;
          size_t table_end = (i == this->total_tasks-1) ? this->n*this->m : (i+1)*task_length;
          size_t row_start = table_start / this->m;
          size_t row_end = (table_end - 1) / this->m;
          this->partical_end_points.emplace_back(std::make_pair(row_start, row_end));
          auto subTask(new TASK<NUMX, NUMY, NUMT, NUMR>(this->optimal_block, i, this->pIdx, this->enc, this->runtime, this->eval));
          subTask->circuit_construct(this->shapeX, this->shapeY, table_start, table_end);
          this->subTasks.emplace_back(subTask);
      }
    }

protected:

    virtual void final_combine(std::vector<std::vector<NUMR>>& partical_result_list, std::vector<std::vector<bool>>& partical_flag_list) override {
        for(int i=0; i<this->total_tasks; i++){
            size_t res_start = this->partical_end_points[i].first, res_end = this->partical_end_points[i].second;
            for(int j=0; j<(res_end - res_start + 1); j++){
                if(partical_flag_list[i][j]) this->res[res_start+j] = partical_result_list[i][j];
                else this->res[res_start+j] = this->res[res_start+j] + partical_result_list[i][j];
            }
        }
        return;
    }
};


int ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc);