#include "CipherIndex.h"
#include "BuildingBlocks.h"

using namespace std;
using namespace oc;
using namespace aby3;

#define SEPARATE_BY_BLOCK_SIZE = true;
static const int BLOCK_SIZE = 5000;
static const int THREAD_NUM = 10;

int cipher_index(u64 pIdx, sf64Matrix<D8> &sharedM, si64Matrix &cipherIndex, sf64Matrix<D8> &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  return 0;
}


int cipher_argsort_offset(int pIdx, si64Matrix& isharedM, si64Matrix& res, Sh3Evaluator& eval, Sh3Runtime& runtime, Sh3Encryptor& enc, Sh3Task& task, int offsetLeft, int offsetRight){
  
  // 1. compute the diff between target elements from offsetLeft to offsetRight.
  u64 block_length = offsetRight - offsetLeft;
  u64 n = isharedM.rows();
  si64Matrix diff(block_length, 1);
  // si64Matrix isharedM = sharedM.i64Cast();

  u64 start_i = offsetLeft / n, start_j = offsetLeft % n;
  u64 end_i = offsetRight / n, end_j = offsetRight % n;
  u64 p = 0, i = start_i, j = start_j;

  while (i < end_i) {
    for (j; j < n; j++) {
      diff.mShares[0](p, 0) = isharedM.mShares[0](j, 0) - isharedM.mShares[0](i, 0);
      diff.mShares[1](p, 0) = isharedM.mShares[1](j, 0) - isharedM.mShares[1](i, 0);
      p++;
    }
    j = 0;
    i++;
  }
  for (j = 0; j < end_j; j++) {
    diff.mShares[0](p, 0) = isharedM.mShares[0](j, 0) - isharedM.mShares[0](end_i, 0);
    diff.mShares[1](p, 0) = isharedM.mShares[1](j, 0) - isharedM.mShares[1](end_i, 0);
    p++;
  }

  // 2. get the msb (comparision result) of the diff computed above.
  sbMatrix pairwise_comp;
  fetch_msb(pIdx, diff, pairwise_comp, eval, runtime, task);

  // 3. convert the comparision results to secret ints.
  i64Matrix ones(block_length, 1);
  for(int i=0; i<block_length; i++) ones(i, 0) = 1;
  si64Matrix ipairwise(block_length, 1);

  // synchronzed initialize the secret ipairwise matrix to ones.
  if(pIdx == 0){
    for(int i=0; i<ipairwise.mShares[0].size(); ++i){
      ipairwise.mShares[0](i) = enc.mShareGen.getShare() + ones(i);
    }
    runtime.mComm.mNext.asyncSendCopy(ipairwise.mShares[0].data(), ipairwise.mShares[0].size());
    auto fu = runtime.mComm.mPrev.asyncRecv(ipairwise.mShares[1].data(), ipairwise.mShares[1].size());
    fu.get();
  }else{
    for(int i=0; i<ipairwise.mShares[0].size(); ++i){
      ipairwise.mShares[0](i) = enc.mShareGen.getShare();
    }
    runtime.mComm.mNext.asyncSendCopy(ipairwise.mShares[0].data(), ipairwise.mShares[0].size());
    auto fu = runtime.mComm.mPrev.asyncRecv(ipairwise.mShares[1].data(), ipairwise.mShares[1].size());
    fu.get();
  }
  // multiply ones with sbMatrix pairwise_comp, convert to secret ints.
  cipher_mul_seq(pIdx, ipairwise, pairwise_comp, ipairwise, eval, enc, runtime);

  // 4. reduce the result to the global res.
  i = start_i, j = start_j;
  p = 0;
  while (i < end_i) {
    for (j; j < n; j++) {
      res.mShares[0](i, 0) += ipairwise.mShares[0](p, 0);
      res.mShares[1](i, 0) += ipairwise.mShares[1](p, 0);
      p++;
    }
    j = 0;
    i++;
  }
  for (j = 0; j < end_j; j++) {
    res.mShares[0](i, 0) += ipairwise.mShares[0](p, 0);
    res.mShares[1](i, 0) += ipairwise.mShares[1](p, 0);
    p++;
  }
  return 0;
}

// repeat-then-reduce version
int rtr_cipher_argsort(int pIdx, si64Matrix& sharedM, si64Matrix& res, Sh3Evaluator& eval, Sh3Runtime&runtime, Sh3Encryptor& enc){

  // 0. define the separate volumn.
  u64 n = sharedM.rows();
  u64 expand_length = n * n;

  // 1. initialize res to all zeros matrix.
  i64Matrix zeros(n, 1);
  for(int i=0; i<n; i++) zeros(i, 0) = 0;
  res.resize(n, 1);
  if(pIdx == 0){
    enc.localIntMatrix(runtime, zeros, res).get();
  }
  else{
    enc.remoteIntMatrix(runtime, res).get();
  }

  // 2. separate the task by the largest volumn.
  int numBlocks = (expand_length + BLOCK_SIZE - 1) / BLOCK_SIZE;
  vector<Sh3Task> vecTasks; 

  for(int curBlock=0; curBlock<numBlocks; curBlock++){
    // compute the offsets.
    int offsetLeft = curBlock * BLOCK_SIZE;
    int offsetRight = (curBlock + 1) * BLOCK_SIZE;
    if (curBlock == numBlocks - 1) {
      offsetRight = (int) expand_length;
    }
    // execute computation of each block in parallel.
    Sh3Task taskCur = runtime.noDependencies();
    string task_id = "task-"+to_string(curBlock)+"-"+to_string(pIdx);
    Sh3Task dep = taskCur.then([&](Sh3Task& self){
      cipher_argsort_offset(pIdx, sharedM, res, eval, runtime, enc, self, offsetLeft, offsetRight);
    }, task_id).getClosure();
    dep.get();
  }

  if(runtime.mIsActive){
    cout << "party - " << pIdx << " INDEED EXISTS ACTIVE TASKS" << endl;
  }
  return 0;
}

int cipher_argsort(int pIdx, si64Matrix &sharedM, si64Matrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Encryptor &enc){
    Sh3Task task = runtime.noDependencies();
    u64 n = sharedM.rows();
    return cipher_argsort_offset(pIdx, sharedM, res, eval, runtime, enc, task, 0, (int)n);
}