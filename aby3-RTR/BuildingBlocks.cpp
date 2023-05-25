#include "BuildingBlocks.h"
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>

#include <aby3/sh3/Sh3BinaryEvaluator.h>
#include <aby3/Circuit/CircuitLibrary.h>
using namespace aby3;
using namespace std;
using namespace oc;

static int BASEPORT=5000;

void distribute_setup(u64 partyIdx, IOService &ios, Sh3Encryptor &enc, Sh3Evaluator &eval,
           Sh3Runtime &runtime) {
  CommPkg comm;
  switch (partyIdx) {
    case 0:
      comm.mNext = Session(ios, "10.0.3.2:1419", SessionMode::Server, "01")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.2:1420", SessionMode::Server, "02")
                       .addChannel();
      break;
    case 1:
      comm.mNext = Session(ios, "10.0.3.4:1421", SessionMode::Server, "12")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.2:1419", SessionMode::Client, "01")
                       .addChannel();
      break;
    default:
      comm.mNext = Session(ios, "10.0.3.2:1420", SessionMode::Client, "02")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.4:1421", SessionMode::Client, "12")
                       .addChannel();
      break;
  }
    // Establishes some shared randomness needed for the later protocols
    enc.init(partyIdx, comm, sysRandomSeed());
    // Establishes some shared randomness needed for the later protocols
    eval.init(partyIdx, comm, sysRandomSeed());
    // Copies the Channels and will use them for later protcols.
    runtime.init(partyIdx, comm);
}

void multi_processor_setup(u64 partyIdx, int rank, IOService &ios, Sh3Encryptor &enc, Sh3Evaluator &eval, Sh3Runtime &runtime){
  CommPkg comm;
  string fport, sport, tport;
  switch (partyIdx) {
    case 0:
      fport = std::to_string(BASEPORT + 3*rank);
      sport = std::to_string(BASEPORT + 3*rank + 1);
      comm.mNext = Session(ios, "10.0.3.8:"+fport, SessionMode::Server, "01")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.8:"+sport, SessionMode::Server, "02")
                       .addChannel();
      break;
    case 1:
      fport = std::to_string(BASEPORT + 3*rank);
      tport = std::to_string(BASEPORT + 3*rank + 2);
      comm.mNext = Session(ios, "10.0.3.4:"+tport, SessionMode::Server, "12")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.8:"+fport, SessionMode::Client, "01")
                       .addChannel();
      break;
    default:
      sport = std::to_string(BASEPORT + 3*rank + 1);
      tport = std::to_string(BASEPORT + 3*rank + 2);
      comm.mNext = Session(ios, "10.0.3.8:"+sport, SessionMode::Client, "02")
                       .addChannel();
      comm.mPrev = Session(ios, "10.0.3.4:"+tport, SessionMode::Client, "12")
                       .addChannel();
      break;
  }
  // Establishes some shared randomness needed for the later protocols
  enc.init(partyIdx, comm, sysRandomSeed());
  // Establishes some shared randomness needed for the later protocols
  eval.init(partyIdx, comm, sysRandomSeed());
  // Copies the Channels and will use them for later protcols.
  runtime.init(partyIdx, comm);
}


void basic_setup(u64 partyIdx, IOService &ios, Sh3Encryptor &enc, Sh3Evaluator &eval,
           Sh3Runtime &runtime) {
  CommPkg comm;
  switch (partyIdx) {
    case 0:
      comm.mNext = Session(ios, "127.0.0.1:1213", SessionMode::Server, "01")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1214", SessionMode::Server, "02")
                       .addChannel();
      break;
    case 1:
      comm.mNext = Session(ios, "127.0.0.1:1215", SessionMode::Server, "12")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1213", SessionMode::Client, "01")
                       .addChannel();
      break;
    default:
      comm.mNext = Session(ios, "127.0.0.1:1214", SessionMode::Client, "02")
                       .addChannel();
      comm.mPrev = Session(ios, "127.0.0.1:1215", SessionMode::Client, "12")
                       .addChannel();
      break;
  }
    // Establishes some shared randomness needed for the later protocols
    enc.init(partyIdx, comm, sysRandomSeed());
    // Establishes some shared randomness needed for the later protocols
    eval.init(partyIdx, comm, sysRandomSeed());
    // Copies the Channels and will use them for later protcols.
    runtime.init(partyIdx, comm);
}


int cipher_mul_seq(int pIdx, const si64Matrix &sharedA, const sbMatrix &sharedB, si64Matrix &res, Sh3Evaluator &eval, Sh3Encryptor& enc, Sh3Runtime &runtime){
  switch (pIdx)
  {
  case 0:
  {
    vector<array<i64, 2>> s0(sharedA.size());
    BitVector c1(sharedA.size());
    for (u64 i = 0; i < s0.size(); ++i)
    {
        auto bb = sharedB.mShares[0](i) ^ sharedB.mShares[1](i);
        auto zeroShare = enc.mShareGen.getShare();

        s0[i][bb] = zeroShare;
        s0[i][bb ^ 1] = sharedA.mShares[1](i) + zeroShare;
        c1[i] = static_cast<u8>(sharedB.mShares[1](i));
    }
    eval.mOtNext.send(runtime.mComm.mNext, s0);
    eval.mOtPrev.send(runtime.mComm.mPrev, s0);

    // share 1: from p1 to p0, p2 
    eval.mOtPrev.help(runtime.mComm.mPrev, c1);
    auto fu1 = runtime.mComm.mPrev.asyncRecv(res.mShares[0].data(), res.size()).share();
    i64* dd = res.mShares[1].data();
    auto fu2 = SharedOT::asyncRecv(runtime.mComm.mNext, runtime.mComm.mPrev, std::move(c1), { dd, i64(res.size()) }).share();
    fu1.get();
    fu2.get();
    break;
  }
  case 1: {
    vector<array<i64, 2>> s1(sharedA.size());
    BitVector c0(sharedA.size());
    for (u64 i = 0; i < s1.size(); ++i)
    {
        auto bb = sharedB.mShares[0](i) ^ sharedB.mShares[1](i);
        auto ss = enc.mShareGen.getShare();

        s1[i][bb] = ss;
        s1[i][bb ^ 1] = (sharedA.mShares[0](i) + sharedA.mShares[1](i)) + ss;
        c0[i] = static_cast<u8>(sharedB.mShares[0](i));
    }
    // share 0: from p0 to p1,p2
    eval.mOtNext.help(runtime.mComm.mNext, c0);

    // share 1: from p1 to p0,p2 
    eval.mOtNext.send(runtime.mComm.mNext, s1);
    eval.mOtPrev.send(runtime.mComm.mPrev, s1);

    // share 0: from p0 to p1,p2
    i64* dd = res.mShares[0].data();
    auto fu1 = SharedOT::asyncRecv(runtime.mComm.mPrev, runtime.mComm.mNext, std::move(c0), { dd, i64(res.size()) }).share();
    // share 1:
    auto fu2 = runtime.mComm.mNext.asyncRecv(res.mShares[1].data(), res.size()).share();
    fu1.get();
    fu2.get();
    break;
  }
  case 2: {
    BitVector c0(sharedA.size()), c1(sharedA.size());
    std::vector<i64> s0(sharedA.size()), s1(sharedA.size());
    for (u64 i = 0; i < sharedA.size(); ++i)
    {
        c0[i] = static_cast<u8>(sharedB.mShares[1](i));
        c1[i] = static_cast<u8>(sharedB.mShares[0](i));
        s0[i] = s1[i] = enc.mShareGen.getShare();
    }
    // share 0: from p0 to p1,p2
    eval.mOtPrev.help(runtime.mComm.mPrev, c0);
    runtime.mComm.mNext.asyncSend(std::move(s0));

    // share 1: from p1 to p0,p2 
    eval.mOtNext.help(runtime.mComm.mNext, c1);
    runtime.mComm.mPrev.asyncSend(std::move(s1));

    // share 0: from p0 to p1,p2
    i64* dd0 = res.mShares[1].data();
    auto fu1 = SharedOT::asyncRecv(runtime.mComm.mNext, runtime.mComm.mPrev, std::move(c0), { dd0, i64(res.size()) }).share();

    // share 1: from p1 to p0,p2
    i64* dd1 = res.mShares[0].data();
    auto fu2 = SharedOT::asyncRecv(runtime.mComm.mPrev, runtime.mComm.mNext, std::move(c1), { dd1, i64(res.size()) }).share();

    fu1.get();
    fu2.get();
    break;
  }
  default:
    throw std::runtime_error(LOCATION);
  }
  return 0;
}


// sequential multiplication between plain i64 and secret bits.
int pi_cb_mul(int pIdx, const i64Matrix &plainA, const sbMatrix &sharedB, si64Matrix &res, Sh3Evaluator &eval, Sh3Encryptor& enc, Sh3Runtime &runtime){

  switch(pIdx) {
    case 0: {
      std::vector<std::array<i64, 2>> s0(plainA.size());
      for(u64 i=0; i<s0.size(); i++){

        auto bb = sharedB.mShares[0](i) ^ sharedB.mShares[1](i);
        auto zeroShare = enc.mShareGen.getShare();

        s0[i][bb] = zeroShare;
        s0[i][bb ^ 1] = plainA(i) + zeroShare;
      }
      eval.mOtNext.send(runtime.mComm.mNext, s0);
      eval.mOtPrev.send(runtime.mComm.mPrev, s0);

      auto fu1 = runtime.mComm.mNext.asyncRecv(res.mShares[0].data(), res.size()).share();
      auto fu2 = runtime.mComm.mPrev.asyncRecv(res.mShares[1].data(), res.size()).share();
      fu1.get();
      fu2.get();
      break;
    }
    case 1:{
      BitVector c0(sharedB.rows());
      for(u64 i=0; i<sharedB.rows(); i++){
        res.mShares[1](i) = enc.mShareGen.getShare();
        c0[i] = static_cast<u8>(sharedB.mShares[0](i));
      }
      // cout << "p1 next no: " << runtime.mComm.mNext << endl;
      eval.mOtNext.help(runtime.mComm.mNext, c0);
      runtime.mComm.mPrev.asyncSendCopy(res.mShares[1].data(), res.size());

      i64* dd = res.mShares[0].data();
      auto fu1 = SharedOT::asyncRecv(runtime.mComm.mPrev, runtime.mComm.mNext, std::move(c0), {dd, i64(res.size())}).share();
      fu1.get();
      break;   
    }
    case 2:{
      BitVector c0(sharedB.rows());
      for(u64 i = 0; i < sharedB.rows(); i++){
        res.mShares[0](i) = enc.mShareGen.getShare();
        c0[i] = static_cast<u8>(sharedB.mShares[1](i));
      }

      // cout << "p2 prev no: " << runtime.mComm.mPrev << endl;
      eval.mOtPrev.help(runtime.mComm.mPrev, c0);
      runtime.mComm.mNext.asyncSendCopy(res.mShares[0].data(), res.size());
      i64* dd0 = res.mShares[1].data();
      auto fu1 = SharedOT::asyncRecv(runtime.mComm.mNext, runtime.mComm.mPrev, std::move(c0), {dd0, i64(res.size())}).share();
      fu1.get();
      break;
    }
    default:
      throw std::runtime_error(LOCATION);
  }

  return 0;
}


// synchronized version of fetch_msb.
int fetch_msb(int pIdx, si64Matrix &diffAB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime, Sh3Task& task){

  // 1. set the binary circuits.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(diffAB.size(), 64);
  circuitInput1.resize(diffAB.size(), 64);

  // 2. let party0 adds up two shares x0 and x1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<diffAB.size(); j++){
        circuitInput0.mShares[0](j) = diffAB.mShares[0](j)+ diffAB.mShares[1](j);
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;  
    }
    case 1: {
      circuitInput1.mShares[0].resize(diffAB.rows(), diffAB.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), diffAB.mShares[0].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(diffAB.rows(), diffAB.cols());
      memcpy(circuitInput1.mShares[1].data(), diffAB.mShares[1].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
    }
  }

  // 3. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  // 4. extract the msb of the difference.
  auto cir = lib.int_comp_helper(sizeof(i64)*8);
  binEng.setCir(cir, diffAB.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);

  // sequentially evaluate the circuit.
  res.resize(diffAB.size(), 1);
  binEng.validateMemory();
  binEng.distributeInputs();
  binEng.roundCallback(runtime.mComm, task);
  binEng.getOutput(0, res);

  return 0;
}


int fetch_msb(int pIdx, si64Matrix &diffAB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // 1. set the binary circuits.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(diffAB.size(), 64);
  circuitInput1.resize(diffAB.size(), 64);

  // 2. let party0 adds up two shares x0 and x1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<diffAB.size(); j++){
        circuitInput0.mShares[0](j) = diffAB.mShares[0](j)+ diffAB.mShares[1](j);
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      circuitInput1.mShares[0].resize(diffAB.rows(), diffAB.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), diffAB.mShares[0].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(diffAB.rows(), diffAB.cols());
      memcpy(circuitInput1.mShares[1].data(), diffAB.mShares[1].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
    }
  }

  // 3. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  // 4. extract the msb of the difference.
  auto cir = lib.int_comp_helper(sizeof(i64)*8);
  binEng.setCir(cir, diffAB.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);

  auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(diffAB.size(), 1);
    binEng.getOutput(0, res);
  });
  dep.get();
  runtime.runUntilTaskCompletes(runtime);

  if(runtime.mIsActive){
    cout << "party - " << pIdx << " INDEED EXISTS ACTIVE TASKS" << endl;
  }

  return 0;
}


int cipher_gt(int pIdx, si64Matrix &sharedA, si64Matrix &sharedB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // comparision between sharedA and sharedB, diffAB = sharedA - sharedB, sub by sint type.
  si64Matrix diffAB = sharedB - sharedA;
  Sh3Task task = runtime.noDependencies();
  return fetch_msb(pIdx, diffAB, res, eval, runtime, task);
}


int cipher_gt(int pIdx, si64Matrix &sharedA, vector<int> &plainB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){
  
  // 1. set the binary circuits.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(sharedA.size(), 64);
  circuitInput1.resize(sharedA.size(), 64);

  // 2. let party0 adds up two shares x0 and x1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<sharedA.size(); j++){
        circuitInput0.mShares[0](j) = plainB[j] - (sharedA.mShares[0](j)+ sharedA.mShares[1](j));
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      circuitInput1.mShares[0].resize(sharedA.rows(), sharedA.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), sharedA.mShares[0].data(), sharedA.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(sharedA.rows(), sharedA.cols());
      memcpy(circuitInput1.mShares[1].data(), sharedA.mShares[1].data(), sharedA.size() * sizeof(i64));

      circuitInput0.mShares[0].setZero();
    }
  }

  // 3. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  // 4. extract the msb of the difference.
  auto cir = lib.int_comp_helper(sizeof(i64)*8);
  binEng.setCir(cir, sharedA.size(), eval.mShareGen);
  binEng.setInput(0, circuitInput0);
  binEng.setInput(1, circuitInput1);

  // sequentially evaluate the circuit.
  Sh3Task task = runtime.noDependencies();
  res.resize(sharedA.size(), 1);
  binEng.validateMemory();
  binEng.distributeInputs();
  binEng.roundCallback(runtime.mComm, task);
  binEng.getOutput(0, res);

  return 0;
}


int cipher_eq(int pIdx, si64Matrix &intA, si64Matrix &intB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // 1. set the difference between AB and BA.
  si64Matrix diffAB = intA - intB;
  si64Matrix diffBA = intB - intA;

  // 2. set the correspondong boolean circuit.
  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;

  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(intA.size(), 64);
  circuitInput1.resize(intB.size(), 64);

  sbMatrix binGt(intA.size(), 1);
  sbMatrix binLt(intB.size(), 1);

  fetch_msb(pIdx, diffAB, binGt, eval, runtime);
  fetch_msb(pIdx, diffBA, binLt, eval, runtime);

  auto cirOr = lib.bits_nor_helper(1);
  binEng.setCir(cirOr, intA.size(), eval.mShareGen);
  binEng.setInput(0, binGt);
  binEng.setInput(1, binLt);
  binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
    res.resize(intA.size(), 1);
    binEng.getOutput(0, res);
  }).get();

  return 0;
}


int circuit_cipher_eq(int pIdx, si64Matrix &intA, si64Matrix &intB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // 1. set the difference between intA and intB
  si64Matrix diffAB = intA - intB;

  // 2. set the correspondong boolean circuit.
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(intA.size(), 64);
  circuitInput1.resize(intB.size(), 64);

  // 3. let party0 adds up two shares x0 + x1 and share with p1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<diffAB.size(); j++){
        circuitInput0.mShares[0](j) = - diffAB.mShares[0](j) - diffAB.mShares[1](j);
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      circuitInput1.mShares[0].resize(diffAB.rows(), diffAB.cols());
      circuitInput1.mShares[1].setZero();
      memcpy(circuitInput1.mShares[0].data(), diffAB.mShares[0].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].resize(diffAB.rows(), diffAB.cols());
      memcpy(circuitInput1.mShares[1].data(), diffAB.mShares[1].data(), diffAB.size() * sizeof(i64));
      circuitInput0.mShares[0].setZero();
    }
  }

  // 4. binary secrets reshare.
  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(), circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(), circuitInput0.mShares[1].size());
  fu.get();

  fetch_eq_res(pIdx, circuitInput0, circuitInput1, res, eval, runtime);

  return 0;
}

int vector_cipher_eq(int pIdx, std::vector<aby3::si64>& intA, std::vector<int>& intB, sbMatrix &res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // set the correspondong boolean circuit.
  sbMatrix circuitInput0;
  sbMatrix circuitInput1;
  circuitInput0.resize(intA.size(), sizeof(aby3::i64) * 8);
  circuitInput1.resize(intB.size(), sizeof(aby3::i64) * 8);

  // 3. let party0 adds up two shares x0 + x1 and share with p1; let party1 and party2 share x2.
  switch(pIdx){
    case 0: {
      for(u64 j=0; j<intA.size(); j++){
        circuitInput0.mShares[0](j) = intB[j] - intA[j].mData[0] - intA[j].mData[1];
      }
      circuitInput1.mShares[0].setZero();
      circuitInput1.mShares[1].setZero();
      break;
    }
    case 1: {
      for(u64 j=0; j<intA.size(); j++){
        circuitInput1.mShares[0](j) = intA[j].mData[0];
      }
      circuitInput1.mShares[1].setZero();
      circuitInput0.mShares[0].setZero();
      break;
    }
    case 2: {
      for(u64 j=0; j<intA.size(); j++){
        circuitInput1.mShares[1](j) = intA[j].mData[1];
      }
      circuitInput1.mShares[0].setZero();
      circuitInput0.mShares[0].setZero();
    }
  }

  runtime.mComm.mNext.asyncSend(circuitInput0.mShares[0].data(),
                                circuitInput0.mShares[0].size());
  auto fu = runtime.mComm.mPrev.asyncRecv(circuitInput0.mShares[1].data(),
                                          circuitInput0.mShares[1].size());
  fu.get();

  fetch_eq_res(pIdx, circuitInput0, circuitInput1, res, eval, runtime);
}


int fetch_eq_res(int pIdx, sbMatrix& circuitA, sbMatrix& circuitB, sbMatrix& res, Sh3Evaluator &eval, Sh3Runtime &runtime){

  // construct the eq circuit.
  int bitSize = circuitA.bitCount();
  int i64Size = circuitA.i64Size();

  Sh3BinaryEvaluator binEng;
  CircuitLibrary lib;

  auto cirEq = lib.int_eq(bitSize);
  binEng.setCir(cirEq, i64Size, eval.mShareGen);
  binEng.setInput(0, circuitA);
  binEng.setInput(1, circuitB);

  // sequentially evaluate the eq circuit.
  Sh3Task task = runtime.noDependencies();
  res.resize(i64Size, 1);
  binEng.validateMemory();
  binEng.distributeInputs();
  binEng.roundCallback(runtime.mComm, task);
  binEng.getOutput(0, res);

  return 0;
}


int init_zeros(int pIdx, Sh3Encryptor &enc, Sh3Runtime &runtime, si64Matrix &res, int n){
    i64Matrix zeros(n, 1);
    for(int i=0; i<n; i++) zeros(i, 0) = 0;
    res.resize(n, 1);
    if(pIdx == 0) {
        enc.localIntMatrix(runtime, zeros, res).get();
    }
    else{
        enc.remoteIntMatrix(runtime, res).get();
    }
    return 0;
}

int init_ones(int pIdx, Sh3Encryptor &enc, Sh3Runtime &runtime, si64Matrix &res, int n){
    i64Matrix ones(n, 1);
    for(int i=0; i<n; i++) ones(i, 0) = 1;
    res.resize(n, 1);
    if(pIdx == 0) {
        enc.localIntMatrix(runtime, ones, res).get();
    }
    else{
        enc.remoteIntMatrix(runtime, res).get();
    }
    return 0;
}
