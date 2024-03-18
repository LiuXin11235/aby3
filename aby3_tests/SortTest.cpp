#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Sort.h"
#include "../aby3-RTR/BuildingBlocks.h"


using namespace oc;
using namespace aby3;

int bc_sort_test(oc::CLP &cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if (role == 0) {
        debug_info("RUN BC Sort TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);
    
    // prepare the data.
    aby3::i64Matrix data(10, 1);
    aby3::i64Matrix res(10, 1);
    for(size_t i=0; i<10; i++) {
        data(i, 0) = 9 - i;
        if(i < 5){
            res(i, 0) = 5 + i;
        }
        else{
            res(i, 0) = i - 5;
        }
    }
    
    aby3::sbMatrix enc_data(10, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(10);
    for(size_t i=0; i<10; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    // sort the data.
    std::vector<size_t> lows = {0, 5};
    std::vector<size_t> highs = {5, 10};

    bc_sort_different(enc_data_vec, lows, highs, role, enc, eval, runtime, 1024);
    aby3::sbMatrix enc_test(10, 64);
    for(size_t i=0; i<10; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    // check the result.
    aby3::i64Matrix test(10, 1);
    enc.revealAll(runtime, enc_test, test).get();


    if(role == 0){
        check_result("BC Sort Test", test, res);
    }

    return 0;
}