#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Network/IOService.h>
#include <mpi.h>

#include "BuildingBlocks.h"
#include "debug.h"
#include "./Pair_then_Reduce/include/datatype.h"
#include "./Pair_then_Reduce/include/tasks.h"

#define OPTIMAL_BLOCK 10000
#define TASKS 20
// #define DEBUG


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

    #ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if(std::is_same<NUMR, aby3::si64>::value){
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "expandX: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_vector(expandX, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "expandY: " << binfo->t_start << std::endl;
      for(int i=0; i<expandY.size(); i++) ofs << expandY[i] << " ";
      ofs << std::endl;
      ofs.close();
      debug_output_vector(expandY, *(this->runtime), *(this->enc));

      aby3::si64Matrix expandXM(block_length);
      for(int i=0; i<block_length; i++) expandXM(i, 0, expandX[i]);
    }
    #endif    

    vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *(this->eval), *(this->runtime));

    #ifdef DEBUG
    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "partTable: " << binfo->t_start << std::endl;
    ofs.close();
    debug_output_matrix(partTable, *(this->runtime), *(this->enc), this->pIdx, *(this->eval));
    #endif


    // pairwise vector abmul.
    aby3::si64Matrix expandV;
    expandV.resize(block_length, 1);
    for(size_t i=0; i<block_length; i++){
      expandV(i, 0, this->selectV[binfo->t_start+i]);
    }
    cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *(this->eval), *(this->enc), *(this->runtime));

    #ifdef DEBUG
    if(std::is_same<NUMR, aby3::si64>::value){
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "table_start: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_matrix(expandV, *(this->runtime), *(this->enc));
    }
    #endif

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

    #ifdef DEBUG
    // debug
    if(std::is_same<NUMR, aby3::si64>::value){
      std::vector<aby3::si64> prob_vec(binfo->res_len);
      for(int i=0; i<binfo->res_len; i++) prob_vec[i] = this->partical_res[i+binfo->r_start];

      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "before add local_res = " << binfo->r_start << " len: " << binfo->res_len << std::endl;
      ofs.close();

      debug_output_vector(prob_vec, *(this->runtime), *(this->enc));
    }
    #endif

    for(int i=0; i<binfo->res_len; i++){
        this->partical_res[i+binfo->r_start] = this->partical_res[i+binfo->r_start] + local_res[i];
    }

    #ifdef DEBUG
    // debug
    if(std::is_same<NUMR, aby3::si64>::value){
      std::vector<aby3::si64> prob_vec(binfo->res_len);
      for(int i=0; i<binfo->res_len; i++) prob_vec[i] = this->partical_res[i+binfo->r_start];

      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "res_start = " << binfo->r_start << " len: " << binfo->res_len << std::endl;
      ofs.close();

      debug_output_vector(prob_vec, *(this->runtime), *(this->enc));
    }
    #endif

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
    // void task_split(size_t table_size, size_t task_length) override {
    //   for(int i=0; i<this->total_tasks; i++){
    //       size_t table_start = i*task_length;
    //       size_t table_end = (i == this->total_tasks-1) ? this->n*this->m : (i+1)*task_length;
    //       size_t row_start = table_start / this->m;
    //       size_t row_end = (table_end - 1) / this->m;
    //       this->partical_end_points.emplace_back(std::make_pair(row_start, row_end));
    //       auto subTask(new TASK<NUMX, NUMY, NUMT, NUMR>(this->optimal_block, i, this->pIdx, this->enc, this->runtime, this->eval));
    //       subTask->circuit_construct(this->shapeX, this->shapeY, table_start, table_end);
    //       for(int j=0; j<subTask->res_length; j++) subTask->partical_res[j] = this->default_value;

    //       if(std::is_same<NUMR, aby3::si64>::value){

    //         std::ofstream ofs(debugFile, std::ios_base::app);
    //         ofs << "test init" << " task_id: " << i << std::endl;
    //         ofs.close();

    //         debug_output_vector(subTask->partical_res, this->runtime, this->enc);
    //       }

    //       this->subTasks.emplace_back(subTask);
    //   }
    // }

    void create_sub_task(size_t optimal_block, int task_id, size_t table_start, size_t table_end) override {
      auto subTask(new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx, this->enc, this->runtime, this->eval));
      subTask->circuit_construct(this->shapeX, this->shapeY, table_start, table_end);
      for(int j=0; j<subTask->res_length; j++) subTask->partical_res[j] = this->default_value;
      this->subTasks.emplace_back(subTask);
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


template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,   
          template<typename, typename, typename, typename> class TASK>
class MPISecretIndex : public MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>{

public:
    // aby3 info
    int pIdx;
    aby3::Sh3Encryptor& enc;
    aby3::Sh3Runtime& runtime;
    aby3::Sh3Evaluator& eval;

    // setup all the aby3 environment variables, pIdx and rank.
    using MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::MPIPTRTask;
    MPISecretIndex(int tasks, size_t optimal_block, const int pIdx, aby3::Sh3Encryptor &enc, aby3::Sh3Runtime& runtime, aby3::Sh3Evaluator& eval):
    pIdx(pIdx), enc(enc), runtime(runtime), eval(eval),
    MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block){
      // std::cout << "in init, i am " << this->rank << std::endl;
    }

    // override sub_task create.
    void create_sub_task(size_t optimal_block, int task_id, size_t table_start, size_t table_end) override {
      auto subTaskPtr(new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx, this->enc, this->runtime, this->eval));
      this->subTask.reset(subTaskPtr);
      this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start, table_end);
      for(int j=0; j<this->subTask->res_length; j++) this->subTask->partical_res[j] = this->default_value;
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
  }
};


int ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc);

int mpi_ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc);