#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Network/IOService.h>
#include <mpi.h>
#include <iomanip>

#include "./Pair_then_Reduce/include/datatype.h"
#include "./Pair_then_Reduce/include/tasks.h"
#include "BuildingBlocks.h"
#include "debug.h"

#define OPTIMAL_BLOCK 100
#define TASKS 5
// #define DEBUG

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class SubIndex : public SubTask<NUMX, NUMY, NUMT, NUMR> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor* enc;
  aby3::Sh3Runtime* runtime;
  aby3::Sh3Evaluator* eval;

  using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

  SubIndex(const size_t optimal_block, const int task_id, const int pIdx,
           aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
           aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
    this->have_selective = true;
  }

  virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "resLeft: " << std::endl;
      ofs.close();
      debug_output_vector(resLeft, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "resRight: " << std::endl;
      ofs.close();
      debug_output_vector(resRight, *(this->runtime), *(this->enc));
    }
#endif
    for (int i = 0; i < resLeft.size(); i++)
      local_res[i] = resLeft[i] + resRight[i];

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "local_res: " << std::endl;
      ofs.close();
      debug_output_vector(local_res, *(this->runtime), *(this->enc));
    }
#endif
    return;
  }

 protected:
  virtual void compute_local_table(std::vector<NUMX>& expandX,
                                   std::vector<NUMY>& expandY,
                                   std::vector<NUMT>& local_table,
                                   BlockInfo* binfo) {
    aby3::u64 block_length = binfo->block_len;
    // pairwise vector equal test.
    aby3::sbMatrix partTable;

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "expandX: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_vector(expandX, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "expandY: " << binfo->t_start << std::endl;
      for (int i = 0; i < expandY.size(); i++) ofs << expandY[i] << " ";
      ofs << std::endl;
      ofs.close();
      // debug_output_vector(expandY, *(this->runtime), *(this->enc));

      // aby3::si64Matrix expandXM(block_length);
      // for(int i=0; i<block_length; i++) expandXM(i, 0, expandX[i]);
    }
#endif

    vector_cipher_eq(this->pIdx, expandX, expandY, partTable, *(this->eval),
                     *(this->runtime));

#ifdef DEBUG
    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "partTable: " << binfo->t_start << std::endl;
    ofs.close();
    debug_output_matrix(partTable, *(this->runtime), *(this->enc), this->pIdx,
                        *(this->eval));
#endif

    // pairwise vector abmul.
    aby3::si64Matrix expandV;
    expandV.resize(block_length, 1);
    for (size_t i = 0; i < block_length; i++) {
      expandV(i, 0, this->selectV[binfo->t_start + i]);
    }
    cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *(this->eval),
                   *(this->enc), *(this->runtime));

#ifdef DEBUG
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "table_start: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_matrix(expandV, *(this->runtime), *(this->enc));
    }
#endif

    // transfor to the local_table
    for (size_t i = 0; i < block_length; i++) local_table[i] = expandV(i, 0);
  }
};

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class SubRank : public SubTask<NUMX, NUMY, NUMT, NUMR> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor* enc;
  aby3::Sh3Runtime* runtime;
  aby3::Sh3Evaluator* eval;

  using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

  SubRank(const size_t optimal_block, const int task_id, const int pIdx,
          aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
          aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
    this->have_selective = false;
  }

  virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "resLeft: " << std::endl;
      ofs.close();
      debug_output_vector(resLeft, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "resRight: " << std::endl;
      ofs.close();
      debug_output_vector(resRight, *(this->runtime), *(this->enc));
    }
#endif
    for (int i = 0; i < resLeft.size(); i++)
      local_res[i] = resLeft[i] + resRight[i];

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "local_res: " << std::endl;
      ofs.close();
      debug_output_vector(local_res, *(this->runtime), *(this->enc));
    }
#endif
    return;
  }

 protected:
  virtual void compute_local_table(std::vector<NUMX>& expandX,
                                   std::vector<NUMY>& expandY,
                                   std::vector<NUMT>& local_table,
                                   BlockInfo* binfo) {
    aby3::u64 block_length = binfo->block_len;

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "expandX: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_vector(expandX, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "expandY: " << binfo->t_start << std::endl;
      for (int i = 0; i < expandY.size(); i++) ofs << expandY[i] << " ";
      ofs << std::endl;
      ofs.close();
      // debug_output_vector(expandY, *(this->runtime), *(this->enc));

      // aby3::si64Matrix expandXM(block_length);
      // for(int i=0; i<block_length; i++) expandXM(i, 0, expandX[i]);
    }
#endif

    vector_cipher_gt(this->pIdx, expandX, expandY, local_table, *(this->eval),
                     *(this->enc), *(this->runtime));

  }
};

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class SubSearch : public SubTask<NUMX, NUMY, NUMT, NUMR> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor* enc;
  aby3::Sh3Runtime* runtime;
  aby3::Sh3Evaluator* eval;

  using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

  SubSearch(const size_t optimal_block, const int task_id, const int pIdx,
          aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
          aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(&enc),
        runtime(&runtime),
        eval(&eval),
        SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
    this->have_selective = true;
  }

  virtual void partical_reduction(std::vector<NUMR>& resLeft,
                                  std::vector<NUMR>& resRight,
                                  std::vector<NUMR>& local_res,
                                  BlockInfo* binfo) override {
#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "resLeft: " << std::endl;
      ofs.close();
      debug_output_vector(resLeft, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "resRight: " << std::endl;
      ofs.close();
      debug_output_vector(resRight, *(this->runtime), *(this->enc));
    }
#endif
    for (int i = 0; i < resLeft.size(); i++)
      local_res[i] = resLeft[i] + resRight[i];

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "local_res: " << std::endl;
      ofs.close();
      debug_output_vector(local_res, *(this->runtime), *(this->enc));
    }
#endif
    return;
  }

 protected:
  virtual void compute_local_table(std::vector<NUMX>& expandX,
                                   std::vector<NUMY>& expandY,
                                   std::vector<NUMT>& local_table,
                                   BlockInfo* binfo) {
    aby3::u64 block_length = binfo->block_len;
    aby3::sbMatrix partTable;

#ifdef DEBUG
    // debug -> aby3 eq has sometimes errorness
    if (std::is_same<NUMR, aby3::si64>::value) {
      std::ofstream ofs(debugFile, std::ios_base::app);
      ofs << "expandX: " << binfo->t_start << std::endl;
      ofs.close();
      debug_output_vector(expandX, *(this->runtime), *(this->enc));

      ofs.open(debugFile, std::ios_base::app);
      ofs << "expandY: " << binfo->t_start << std::endl;
      for (int i = 0; i < expandY.size(); i++) ofs << expandY[i] << " ";
      ofs << std::endl;
      ofs.close();
      // debug_output_vector(expandY, *(this->runtime), *(this->enc));
      // aby3::si64Matrix expandXM(block_length);
      // for(int i=0; i<block_length; i++) expandXM(i, 0, expandX[i]);
    }
#endif

    vector_cipher_ge(this->pIdx, expandX, expandY, partTable, *(this->eval),
                     *(this->enc), *(this->runtime));
    
    // pairwise vector abmul.
    aby3::si64Matrix expandV;
    expandV.resize(block_length, 1);
    for (size_t i = 0; i < block_length; i++) expandV(i, 0, this->selectV[binfo->t_start + i]);

    cipher_mul_seq(this->pIdx, expandV, partTable, expandV, *(this->eval),
                   *(this->enc), *(this->runtime));
    
    for (size_t i = 0; i < block_length; i++) local_table[i] = expandV(i, 0);
  }
};

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class SecretIndex : public PTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  using PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::PTRTask;
  SecretIndex(int tasks, size_t optimal_block, const int pIdx,
              aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
              aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  SecretIndex(const int pIdx, aby3::Sh3Encryptor& enc,
              aby3::Sh3Runtime& runtime, aby3::Sh3Evaluator& eval,
              std::string deployment_profile_name = "./Config/profile.json")
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(deployment_profile_name) {}

  using PTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::task_split;
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTask(new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id,
                                                  this->pIdx, this->enc,
                                                  this->runtime, this->eval));
    subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                               table_end);
    subTask->initial_value = this->default_value;
    for (int j = 0; j < this->n; j++) subTask->res[j] = this->default_value;
    this->subTasks.emplace_back(subTask);
  }

  // void circuit_evaluate(NUMX* dataX, NUMY* dataY, NUMR* selectV, NUMR* res){

  //     // prepare data structures.
  //     this->inputX = fake_repeat(dataX, this->shapeX, this->m, 0);
  //     this->inputY = fake_repeat(dataY, this->shapeY, this->n, 1);
  //     this->res = res;

  //     // compute functions => currently, using sequential.
  //     for(int i=0; i<this->total_tasks; i++){
  //         // call the corresponding functions on different machines.
  //         this->subTasks[i]->circuit_evaluate(this->inputX, this->inputY,
  //         this->selectV);
  //     }

  //     #ifdef DEBUG
  //     // debug -> aby3 eq has sometimes errorness
  //     for(int i=0; i<this->total_tasks; i++){
  //       if(std::is_same<NUMR, aby3::si64>::value){
  //         std::ofstream ofs(debugFile, std::ios_base::app);
  //         ofs << "subTask-" << i << "res: " << std::endl;
  //         ofs.close();
  //         debug_output_vector(this->subTasks[i]->res, this->runtime,
  //         this->enc);

  //         // ofs.open(debugFile, std::ios_base::app);
  //         // ofs << "expandY: " << binfo->t_start << std::endl;
  //         // for(int i=0; i<expandY.size(); i++) ofs << expandY[i] << " ";
  //         // ofs << std::endl;
  //         // ofs.close();
  //         // debug_output_vector(expandY, *(this->runtime), *(this->enc));

  //         // aby3::si64Matrix expandXM(block_length);
  //         // for(int i=0; i<block_length; i++) expandXM(i, 0, expandX[i]);
  //       }
  //     }
  //     #endif

  //     // simulate
  //     for(int i=this->total_tasks-1; i>=0; i--){
  //         size_t left_tasks = this->total_tasks;
  //         size_t send_start = (left_tasks + 1) / 2;
  //         while(i < send_start && left_tasks > 1){
  //             size_t receive_target = i + send_start;
  //             if(receive_target < left_tasks){
  //                 std::vector<NUMR> receive_res =
  //                 this->subTasks[receive_target]->res;
  //                 this->subTasks[i]->partical_reduction(receive_res,
  //                 this->subTasks[i]->res, this->subTasks[i]->res, nullptr);
  //             }

  //             left_tasks = send_start;
  //             send_start = (left_tasks + 1) / 2;
  //         }
  //         if(i >= send_start) size_t end_target = i - send_start;
  //     }
  //     std::copy(this->subTasks[0]->res.begin(), this->subTasks[0]->res.end(),
  //     res);
  // }
};

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class MPISecretIndex : public MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  // setup all the aby3 environment variables, pIdx and rank.
  using MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::MPIPTRTask;
  MPISecretIndex(int tasks, size_t optimal_block, const int pIdx,
                 aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
                 aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  // override sub_task create.
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTaskPtr(
        new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx,
                                         this->enc, this->runtime, this->eval));
    this->subTask.reset(subTaskPtr);
    this->subTask->initial_value = this->default_value;
    this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                                     table_end);
    for (int j = 0; j < this->n; j++)
      this->subTask->res[j] = this->default_value;
  }
};

template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class MPIRank : public MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  // setup all the aby3 environment variables, pIdx and rank.
  using MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::MPIPTRTask;
  MPIRank(int tasks, size_t optimal_block, const int pIdx,
                 aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
                 aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  // override sub_task create.
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTaskPtr(
        new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx,
                                         this->enc, this->runtime, this->eval));
    this->subTask.reset(subTaskPtr);
    this->subTask->initial_value = this->default_value;
    this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                                     table_end);
    for (int j = 0; j < this->n; j++)
      this->subTask->res[j] = this->default_value;
  }
};


template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class MPISearch : public MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  // setup all the aby3 environment variables, pIdx and rank.
  using MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::MPIPTRTask;
  MPISearch(int tasks, size_t optimal_block, const int pIdx,
                 aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
                 aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  // override sub_task create.
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTaskPtr(
        new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx,
                                         this->enc, this->runtime, this->eval));
    this->subTask.reset(subTaskPtr);
    this->subTask->initial_value = this->default_value;
    this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                                     table_end);
    for (int j = 0; j < this->n; j++)
      this->subTask->res[j] = this->default_value;
  }
};


int ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                     std::vector<aby3::si64>& secretIndex,
                     std::vector<aby3::si64>& res, aby3::Sh3Evaluator& eval,
                     aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor& enc);

int mpi_ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                         std::vector<aby3::si64>& secretIndex,
                         std::vector<aby3::si64>& res, aby3::Sh3Evaluator& eval,
                         aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor& enc,
                         int task_num, int opt_B);
