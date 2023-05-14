#include "PTRFunction.h"

int ptr_secret_index(int pIdx, std::vector<aby3::si64>& sharedM,
                 std::vector<aby3::si64>& secretIndex, std::vector<aby3::si64>& res,
                 aby3::Sh3Evaluator &eval, aby3::Sh3Runtime &runtime,
                 aby3::Sh3Encryptor &enc){
  auto ptrTask = new SecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(TASKS, OPTIMAL_BLOCK, pIdx, enc, runtime, eval);
  size_t n = secretIndex.size(), m = sharedM.size();
  int* range_index = new int[m];
  for(int i=0; i<m; i++) range_index[i] = i;

  ptrTask->circuit_construct({n}, {m});
  ptrTask->set_selective_value(sharedM.data(), 0);
  ptrTask->circuit_evaluate(secretIndex.data(), range_index, sharedM.data(), res.data());
}