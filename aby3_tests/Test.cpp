#include "Test.h"
#include <chrono>
#include <thread>
#include <random>
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-RTR/debug.h"
#include "../aby3-Basic/Basics.h"

using namespace oc;
using namespace aby3;
using namespace std;

const int TEST_SIZE = 5;

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