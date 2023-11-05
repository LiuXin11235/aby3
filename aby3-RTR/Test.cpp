#include "Test.h"
#include <chrono>
#include <thread>
#include <random>
#include "BuildingBlocks.h"
#include "debug.h"

using namespace oc;
using namespace aby3;
using namespace std;

const int TEST_SIZE = 50;

bool check_result(const std::string& func_name, i64Matrix& test, i64Matrix& res){
    int size = test.rows();
    bool check_flag = true;
    for(int i=0; i<size; i++){
        if(test(i, 0) != res(i, 0)) check_flag = false;
    }

    if(!check_flag){
        debug_info(func_name + " ERROR !");
        debug_info("test: ");
        debug_output_matrix(test);
        debug_info("result: ");
        debug_output_matrix(res);
    }

    return check_flag;
}

int basic_test(CLP& cmd){

    // get the configs.
    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    distribute_setup((u64)role, ios, enc, eval, runtime);

    // generate test data.
    i64Matrix input_x(TEST_SIZE, 1);
    i64Matrix input_y(TEST_SIZE, 1);
    i64Matrix res_gt(TEST_SIZE, 1);
    i64Matrix res_ge(TEST_SIZE, 1);
    i64Matrix res_eq(TEST_SIZE, 1);
    i64Matrix res_mul(TEST_SIZE, 1);

    for(int i=0; i<TEST_SIZE; i++){
        input_x(i, 0) = i;
        input_y(i, 0) = TEST_SIZE - i;

        res_gt(i, 0) = i > (TEST_SIZE - i);
        res_ge(i, 0) = i >= (TEST_SIZE - i);
        res_mul(i, 0) = i * (TEST_SIZE - i);
        res_eq(i, 0) = i == (TEST_SIZE - i);
    }

    // encrypt the inputs.
    si64Matrix sharedX(TEST_SIZE, 1);
    si64Matrix sharedY(TEST_SIZE, 1);
    if (role == 0) {
        enc.localIntMatrix(runtime, input_x, sharedX).get();
        enc.localIntMatrix(runtime, input_y, sharedY).get();
    } else {
        enc.remoteIntMatrix(runtime, sharedX).get();
        enc.remoteIntMatrix(runtime, sharedY).get();
    }

    // compute for the result.
    sbMatrix shared_gt(TEST_SIZE, 1);
    sbMatrix shared_ge(TEST_SIZE, 1);
    sbMatrix shared_eq(TEST_SIZE, 1);
    si64Matrix shared_mul(TEST_SIZE, 1);

    // if(role == 0){
    //     debug_info("inputs");
    //     debug_info("sharedX: ");
    //     debug_output_matrix(sharedX, runtime, enc);
    //     debug_info("sharedY: ");
    //     debug_output_matrix(sharedY, runtime, enc);
    // }
    i64Matrix plain_x(TEST_SIZE, 1);
    i64Matrix plain_y(TEST_SIZE, 1);
    enc.revealAll(runtime, sharedX, plain_x).get();
    enc.revealAll(runtime, sharedY, plain_y).get();
    if(role == 0){
        debug_info("inputs");
        debug_info("sharedX: ");
        debug_output_matrix(plain_x);
        debug_info("sharedY: ");
        debug_output_matrix(plain_y);
    }

    eval.asyncMul(runtime, sharedY, sharedX, shared_mul).get();
    cipher_gt(role, sharedX, sharedY, shared_gt, eval, runtime);
    circuit_cipher_eq(role, sharedX, sharedY, shared_eq, eval, runtime);
    // cipher_mul_seq(role, sharedX, sharedY, shared_mul, eval, enc, runtime);
    cipher_ge(role, sharedX, sharedY, shared_eq, eval, enc, runtime);

    // check the result.
    i64Matrix test_gt(TEST_SIZE, 1); i64Matrix test_ge(TEST_SIZE, 1);
    i64Matrix test_eq(TEST_SIZE, 1); i64Matrix test_mul(TEST_SIZE, 1);

    enc.revealAll(runtime, shared_gt, test_gt).get();
    enc.revealAll(runtime, shared_ge, test_ge).get();
    enc.revealAll(runtime, shared_eq, test_eq).get();
    enc.revealAll(runtime, shared_mul, test_mul).get();

    if(role == 0){
        check_result("gt", test_gt, res_gt);
        check_result("ge", test_ge, res_ge);
        check_result("eq", test_eq, res_eq);
        check_result("mul", test_mul, res_mul);
    }
    return 0;
}