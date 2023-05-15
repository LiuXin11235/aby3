#include "PTRFunction.h"

int ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc){
  auto ptrTask = new SecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(TASKS, OPTIMAL_BLOCK, pIdx, enc, runtime, eval);
  size_t n = secretIndex.size(), m = sharedM.size();
  int* range_index = new int[m];
  for(int i=0; i<m; i++) range_index[i] = i;
  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  // if(pIdx == 0) enc.localInt(runtime, (aby3::i64)0, dval).get();
  // else enc.remoteInt(runtime, dval).get();
  ptrTask->set_default_value(dval);
  ptrTask->circuit_construct({n}, {m});
  ptrTask->set_selective_value(sharedM.data(), 0);
  ptrTask->circuit_evaluate(secretIndex.data(), range_index, sharedM.data(), res.data());
}

int mpi_ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc){

  auto mpiPtrTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(TASKS, OPTIMAL_BLOCK, pIdx, enc, runtime, eval);

  size_t n = secretIndex.size(), m = sharedM.size();
  int* range_index = new int[m];
  for(int i=0; i<m; i++) range_index[i] = i;
  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;

  // distributed task.
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({n}, {m});
  mpiPtrTask->set_selective_value(sharedM.data(), 0);
  mpiPtrTask->circuit_evaluate(secretIndex.data(), range_index, sharedM.data(), res.data());
  return 0;
}