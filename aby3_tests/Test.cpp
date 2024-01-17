#include "Test.h"
#include <chrono>
#include <thread>
#include <random>
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"
#include "../aby3-Basic/Basics.h"
#include "../aby3-Basic/Shuffle.h"

using namespace oc;
using namespace aby3;
using namespace std;

const int TEST_SIZE = 100;
const int TEST_UNIT_SIZE = 10;

bool check_result(const std::string& func_name, i64Matrix& test, i64Matrix& res){
    int size = test.rows();
    bool check_flag = true;
    for(int i=0; i<size; i++){
        if(test(i, 0) != res(i, 0)) check_flag = false;
    }

    if(!check_flag){
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }

    return check_flag;
}


int bool_basic_test(CLP& cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if(role == 0) {
        debug_info("RUN BOOL TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate test data.
    i64Matrix input_x(TEST_SIZE, 1);
    i64Matrix input_y(TEST_SIZE, 1);
    i64Matrix input_z(TEST_SIZE, 1);
    i64Matrix res_xor(TEST_SIZE, 1);
    i64Matrix res_gt(TEST_SIZE, 1);
    i64Matrix res_eq(TEST_SIZE, 1);
    i64Matrix res_add(TEST_SIZE, 1);    
    i64Matrix res_or(TEST_SIZE, 1);
    i64Matrix res_and(TEST_SIZE, 1);

    for(int i=0; i<TEST_SIZE; i++){
        input_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;
        input_z(i, 0) = -i;
        res_xor(i, 0) = input_x(i, 0) ^ input_y(i, 0);
        res_or(i, 0) = input_x(i, 0) | input_y(i, 0);
        res_and(i, 0) = input_x(i, 0) & input_y(i, 0);
        res_gt(i, 0) = i > (TEST_SIZE - i);
        res_eq(i, 0) = i == (TEST_SIZE - i);
        res_add(i, 0) = TEST_SIZE;
    }

    // encrypt the inputs.
    sbMatrix bsharedX(TEST_SIZE, 1);    
    sbMatrix bsharedY(TEST_SIZE, 1);
    sbMatrix bsharedZ(TEST_SIZE, 1);
    bsharedX.resize(TEST_SIZE, 64);
    bsharedY.resize(TEST_SIZE, 64);

    if (role == 0) {
        enc.localBinMatrix(runtime, input_x, bsharedX).get();
        enc.localBinMatrix(runtime, input_y, bsharedY).get();
        enc.localBinMatrix(runtime, input_z, bsharedZ).get();
    } else {
        enc.remoteBinMatrix(runtime, bsharedX).get();
        enc.remoteBinMatrix(runtime, bsharedY).get();
        enc.remoteBinMatrix(runtime, bsharedZ).get();
    }

    // std::ofstream ofs(debugFile, std::ios_base::app);
    // ofs << "z: " << std::endl;
    // for(int i=0; i<input_z.rows(); i++){
    //     printBits(input_z(i), ofs);
    // }
    // ofs << "x: " << std::endl;
    // for(int i=0; i<input_z.rows(); i++){
    //     printBits(input_x(i), ofs);
    // }

    // xor
    sbMatrix shared_xor(TEST_SIZE, 1);
    i64Matrix test_xor(TEST_SIZE, 1);
    for(int i=0; i<TEST_SIZE; i++){
        shared_xor.mShares[0](i) = bsharedX.mShares[0](i) ^ bsharedY.mShares[0](i);
        shared_xor.mShares[1](i) = bsharedX.mShares[1](i) ^ bsharedY.mShares[1](i);
    }
    enc.revealAll(runtime, shared_xor, test_xor).get();

    // or
    sbMatrix shared_or(TEST_SIZE, 1);
    i64Matrix test_or(TEST_SIZE, 1);
    bool_cipher_or(role, bsharedX, bsharedY, shared_or, enc, eval, runtime);
    enc.revealAll(runtime, shared_or, test_or).get();

    if(role == 0){
        debug_output_matrix(test_or);
        debug_output_matrix(res_or);
    }

    // and
    sbMatrix shared_and(TEST_SIZE, 1);
    i64Matrix test_and(TEST_SIZE, 1);
    bool_cipher_and(role, bsharedX, bsharedY, shared_and, enc, eval, runtime);
    enc.revealAll(runtime, shared_and, test_and).get();

    // gt
    sbMatrix shared_gt(TEST_SIZE, 1);
    i64Matrix test_gt(TEST_SIZE, 1);
    bool_cipher_lt(role, bsharedY, bsharedX, shared_gt, enc, eval, runtime);
    enc.revealAll(runtime, shared_gt, test_gt).get();

    // eq 
    sbMatrix shared_eq(TEST_SIZE, 1);
    i64Matrix test_eq(TEST_SIZE, 1);
    bool_cipher_eq(role, bsharedY, bsharedX, shared_eq, enc, eval, runtime);
    enc.revealAll(runtime, shared_eq, test_eq).get();

    // add
    sbMatrix shared_add(TEST_SIZE, 1);
    i64Matrix test_add(TEST_SIZE, 1);
    bool_cipher_add(role, bsharedX, bsharedY, shared_add, enc, eval, runtime);
    enc.revealAll(runtime, shared_add, test_add).get();

    // check the result.
    if(role == 0){
        check_result("xor", test_xor, res_xor);
        check_result("or", test_or, res_or);
        check_result("gt", test_gt, res_gt);
        check_result("eq", test_eq, res_eq);
        check_result("add", test_add, res_add);
        check_result("and", test_and, res_and);
    }

    return 0;
}


int arith_basic_test(CLP& cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if(role == 0) {
        debug_info("RUN ARITH TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    // distribute_setup((u64)role, ios, enc, eval, runtime);
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate test data.
    i64Matrix input_x(TEST_SIZE, 1);
    i64Matrix input_y(TEST_SIZE, 1);
    f64Matrix<D8> finput_x(TEST_SIZE, 1);
    i64Matrix res_gt(TEST_SIZE, 1);
    i64Matrix res_ge(TEST_SIZE, 1);
    i64Matrix res_eq(TEST_SIZE, 1);
    i64Matrix res_mul(TEST_SIZE, 1);
    i64Matrix res_ab_mul(TEST_SIZE, 1);
    i64Matrix res_ib_mul(TEST_SIZE, 1);
    f64Matrix<D8> res_f_mul(TEST_SIZE, 1);
    f64Matrix<D8> res_fb_mul(TEST_SIZE, 1);


    for(int i=0; i<TEST_SIZE; i++){
        input_x(i, 0) = i;
        finput_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;

        res_gt(i, 0) = i > (TEST_SIZE - i);
        res_ge(i, 0) = i >= (TEST_SIZE - i);
        res_mul(i, 0) = i * (TEST_SIZE - i);
        res_eq(i, 0) = i == (TEST_SIZE - i);
        res_ab_mul(i, 0) = res_gt(i, 0) * input_x(i, 0);
        res_ib_mul(i, 0) = i * res_gt(i, 0);
        res_f_mul(i, 0) = i * i;
        res_fb_mul(i, 0) = i * res_gt(i, 0);
    }

    // encrypt the inputs.
    si64Matrix sharedX(TEST_SIZE, 1);
    si64Matrix sharedY(TEST_SIZE, 1);
    sf64Matrix<D8> fsharedX(TEST_SIZE, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_x, sharedX).get();
        enc.localIntMatrix(runtime, input_y, sharedY).get();
        enc.localFixedMatrix(runtime, finput_x, fsharedX).get();
    } else {
        enc.remoteIntMatrix(runtime, sharedX).get();
        enc.remoteIntMatrix(runtime, sharedY).get();
        enc.remoteFixedMatrix(runtime, fsharedX).get();
    }

    // compute for the result.
    sbMatrix shared_gt(TEST_SIZE, 1);
    sbMatrix shared_ge(TEST_SIZE, 1);
    sbMatrix shared_eq(TEST_SIZE, 1);

    si64Matrix shared_mul(TEST_SIZE, 1);
    si64Matrix shared_ab_mul(TEST_SIZE, 1);
    si64Matrix shared_ib_mul(TEST_SIZE, 1);
    sf64Matrix<D8> shared_f_mul(TEST_SIZE, 1);
    sf64Matrix<D8> shared_fb_mul(TEST_SIZE, 1);

    cipher_gt(role, sharedX, sharedY, shared_gt, eval, runtime);
    circuit_cipher_eq(role, sharedX, sharedY, shared_eq, eval, runtime);
    cipher_ge(role, sharedX, sharedY, shared_ge, eval, enc, runtime);

    cipher_mul(role, sharedX, sharedY, shared_mul, eval, enc, runtime);
    cipher_mul(role, sharedX, shared_gt, shared_ab_mul, eval, enc, runtime);
    cipher_mul(role, input_x, shared_gt, shared_ib_mul, eval, enc, runtime);
    // std::cout << "mul3" << std::endl;
    cipher_mul<D8>(role, fsharedX, fsharedX, shared_f_mul, eval, enc, runtime);
    cipher_mul<D8>(role, fsharedX, shared_gt, shared_fb_mul, eval, enc, runtime);

    // check the result.
    i64Matrix test_gt(TEST_SIZE, 1); i64Matrix test_ge(TEST_SIZE, 1);   
    i64Matrix test_eq(TEST_SIZE, 1); i64Matrix test_mul(TEST_SIZE, 1);
    i64Matrix test_ab_mul(TEST_SIZE, 1); i64Matrix test_ib_mul(TEST_SIZE, 1);
    f64Matrix<D8> test_f_mul(TEST_SIZE, 1); f64Matrix<D8> test_fb_mul(TEST_SIZE, 1); 


    enc.revealAll(runtime, shared_gt, test_gt).get();
    enc.revealAll(runtime, shared_ge, test_ge).get();
    enc.revealAll(runtime, shared_eq, test_eq).get();
    enc.revealAll(runtime, shared_mul, test_mul).get();
    enc.revealAll(runtime, shared_ab_mul, test_ab_mul).get();
    enc.revealAll(runtime, shared_ib_mul, test_ib_mul).get();
    enc.revealAll(runtime, shared_f_mul, test_f_mul).get();
    enc.revealAll(runtime, shared_fb_mul, test_fb_mul).get();


    if(role == 0){
        check_result("gt", test_gt, res_gt);
        check_result("ge", test_ge, res_ge);
        check_result("eq", test_eq, res_eq);
        check_result("mul", test_mul, res_mul);
        check_result("mul_ab", test_ab_mul, res_ab_mul);
        check_result("mul_ib", test_ib_mul, res_ib_mul);
        check_result("mul_f", test_f_mul, res_f_mul);
        check_result("mul_fb", test_fb_mul, res_fb_mul);
    }
    return 0;
}


int initialization_test(CLP& cmd){
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if(role == 0) {
        debug_info("RUN INIT TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // 1. test the correlated randomness.
    oc::block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    oc::block prevSeed = enc.mShareGen.mPrevCommon.getSeed();

    // send the nextSeed to the nextParty.  
    runtime.mComm.mNext.asyncSendCopy(nextSeed);

    // send the prevSeed to the prevParty.
    runtime.mComm.mPrev.asyncSendCopy(prevSeed);

    // get the nextParty's prev seed.
    oc::block nextP_seed;
    runtime.mComm.mNext.recv(nextP_seed);

    // get the prevParty's next seed.
    oc::block prevP_seed;
    runtime.mComm.mPrev.recv(prevP_seed);

    // if(role == 0){
    //     debug_info("next seed: ");
    //     debug_info(nextSeed);
    //     debug_info(nextP_seed);
    //     debug_info("prev seed: ");
    //     debug_info(prevSeed);
    //     debug_info(prevP_seed);
    // }

    // check the seeds.
    if(nextSeed != nextP_seed){
        debug_info("\033[31m P" + to_string(role) + " check: nextSeed ERROR!" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m P" + to_string(role) + " check: nextSeed SUCCESS!" + "\033[0m\n");
    }

    if(prevSeed != prevP_seed){
        debug_info("\033[31m P" + to_string(role) + " check: prevSeed ERROR!" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m P" + to_string(role) + " check: prevSeed SUCCESS!" + "\033[0m\n");
    }

    return 0;
}


int shuffle_test(oc::CLP& cmd){
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if(role == 0) {
        debug_info("RUN SHUFFLE TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // generate the test data.
    std::vector<i64Matrix> input_x(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        input_x[i].resize(TEST_UNIT_SIZE, 1);
        for(size_t j=0; j<TEST_UNIT_SIZE; j++){
            input_x[i](j, 0) = i;
        }
    }

    // encrypt the inputs.
    std::vector<sbMatrix> bsharedX(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        bsharedX[i].resize(TEST_UNIT_SIZE, 1);
        if(role == 0){
            enc.localBinMatrix(runtime, input_x[i], bsharedX[i]).get();
        }
        else{
            enc.remoteBinMatrix(runtime, bsharedX[i]).get();
        }
    }

    // generate the permutation.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = TEST_SIZE;
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    std::vector<size_t> other_permutation(len);
    runtime.mComm.mPrev.asyncSendCopy(next_permutation.data(), next_permutation.size());
    runtime.mComm.mNext.recv(other_permutation.data(), other_permutation.size());

    std::vector<size_t> final_permutation;  
    std::vector<std::vector<size_t>> permutation_list;
    if(role == 0){
        permutation_list = {next_permutation, prev_permutation, other_permutation};
    }
    if(role == 1){
        permutation_list = {prev_permutation, other_permutation, next_permutation};
    }
    if(role == 2){
        permutation_list = {other_permutation, next_permutation, prev_permutation};
    }
    combine_permutation(permutation_list, final_permutation);

    std::vector<i64Matrix> shuffle_res(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        shuffle_res[i].resize(TEST_UNIT_SIZE, 1);
        shuffle_res[i] = input_x[i];
    }
    plain_permutate(final_permutation, shuffle_res);

    std::vector<sbMatrix> bsharedShuffle(TEST_SIZE);

    efficient_shuffle(bsharedX, role, bsharedShuffle, enc, eval, runtime);

    std::vector<i64Matrix> test_res(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        test_res[i].resize(TEST_UNIT_SIZE, 1);
        enc.revealAll(runtime, bsharedShuffle[i], test_res[i]).get();
    }

    // check the shuffle result.
    if(role == 0){
        bool check_flag = true;
        for(size_t i=0; i<TEST_SIZE; i++){
            for(size_t j=0; j<TEST_UNIT_SIZE; j++){
                if(test_res[i](j, 0) != shuffle_res[i](j, 0)){
                    check_flag = false;
                }
            }
        }
        if(check_flag){
            debug_info("\033[32m SHUFFLE CHECK SUCCESS ! \033[0m\n");
        }
        else{
            debug_info("\033[31m SHUFFLE CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for(size_t i=0; i<TEST_SIZE; i++){
                debug_output_matrix(shuffle_res[i]);
            }
            debug_info("Func result: \n");
            for(size_t i=0; i<TEST_SIZE; i++){
                debug_output_matrix(test_res[i]);
            }
        }
    }


    // shuffle with permutation test.
    std::vector<si64> shared_permutation(TEST_SIZE);
    efficient_shuffle_with_random_permutation(bsharedX, role, bsharedShuffle, shared_permutation, enc, eval, runtime);

    std::vector<i64Matrix> test_res2(TEST_SIZE);
    for(size_t i=0; i<TEST_SIZE; i++){
        test_res2[i].resize(TEST_UNIT_SIZE, 1);
        enc.revealAll(runtime, bsharedShuffle[i], test_res2[i]).get();
    }
    sbMatrix shared_permutation_matrix(TEST_SIZE, 1);
    for(size_t i=0; i<TEST_SIZE; i++){
        shared_permutation_matrix.mShares[0](i) = shared_permutation[i].mData[0];
        shared_permutation_matrix.mShares[1](i) = shared_permutation[i].mData[1];
    }
    i64Matrix test_permutation(TEST_SIZE, 1);
    enc.revealAll(runtime, shared_permutation_matrix, test_permutation).get();

    bool check_shuffle = true;
    bool check_permutation = true;

    for(size_t i=0; i<TEST_SIZE; i++){
        if(test_permutation(i, 0) != final_permutation[i]){
            check_permutation = false;
        }
        for(size_t j=0; j<TEST_UNIT_SIZE; j++){
            if(test_res2[i](j, 0) != shuffle_res[i](j, 0)){
                check_shuffle = false;
            }
        }
    }

    // check the final result.
    if(role == 0){
        if(check_shuffle){
            debug_info("\033[32m SHUFFLE in shuffle and permutation CHECK SUCCESS ! \033[0m\n");
        }
        else{
            debug_info("\033[31m SHUFFLE in in shuffle and permutation CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            for(size_t i=0; i<TEST_SIZE; i++){
                debug_output_matrix(shuffle_res[i]);
            }
            debug_info("Func result: \n");
            for(size_t i=0; i<TEST_SIZE; i++){
                debug_output_matrix(test_res2[i]);
            }
        }

        if(check_permutation){
            debug_info("\033[32m PERMUTATION in shuffle and permutation CHECK SUCCESS ! \033[0m\n");
        }
        else{
            debug_info("\033[31m PERMUTATION in in shuffle and permutation CHECK ERROR ! \033[0m\n");
            debug_info("True result: \n");
            debug_output_vector(final_permutation);

            debug_info("Func result: \n");
            debug_output_matrix(test_permutation);
        }
    }
    return 0;
}


int correlation_test(oc::CLP& cmd){
    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    if(role == 0) {
        debug_info("RUN CORRELATION TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    i64Matrix prevMask(TEST_SIZE, 1);
    i64Matrix nextMask(TEST_SIZE, 1);

    // generate the permutation.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();

    get_random_mask(role, prevMask, prevSeed);
    get_random_mask(role, nextMask, nextSeed);

    // check.
    i64Matrix nextP_mask(TEST_SIZE, 1);
    i64Matrix prevP_mask(TEST_SIZE, 1);

    runtime.mComm.mNext.asyncSendCopy(nextMask.data(), nextMask.size());
    runtime.mComm.mPrev.asyncSendCopy(prevMask.data(), prevMask.size());

    runtime.mComm.mNext.recv(nextP_mask.data(), nextP_mask.size());
    runtime.mComm.mPrev.recv(prevP_mask.data(), prevP_mask.size());

    bool next_check = true;
    bool prev_check = true;

        
    for(size_t i=0; i<TEST_SIZE; i++){
        if(nextMask(i, 0) != nextP_mask(i, 0)){
            debug_info(to_string(nextMask(i, 0)) + " " + to_string(nextP_mask(i, 0)) + "\n");
            next_check = false;
        }
        if(prevMask(i, 0) != prevP_mask(i, 0)){
            debug_info(to_string(prevMask(i, 0)) + " " + to_string(prevP_mask(i, 0)) + "\n");
            prev_check = false;
        }
    }

    if(!next_check){
        debug_info("\033[31m P" + to_string(role) + " check: next randomness ERROR!" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m P" + to_string(role) + " check: next randomness SUCCESS!" + "\033[0m\n");
    }

    if(!prev_check){
        debug_info("\033[31m P" + to_string(role) + " check: prev randomness ERROR!" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m P" + to_string(role) + " check: prev randomness SUCCESS!" + "\033[0m\n");
    }

    return 0;
}