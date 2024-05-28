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
        bool check_flag = check_result("BC Sort Test", test, res);
        if(!check_flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }

    return 0;
}

int bc_sort_corner_test(oc::CLP &cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN BC Sort Corner TEST");
    }

    // prepare the data.
    aby3::i64Matrix data(10, 1);
    std::vector<int> res_vec = {9, 7, 8, 5, 6, 2, 3, 4, 0, 1};
    aby3::i64Matrix res(10, 1);
    for(size_t i=0; i<10; i++) {
        data(i, 0) = 9 - i;
        res(i, 0) = res_vec[i];
    }

    // sort the data.
    std::vector<size_t> lows = {1, 3, 5, 8};
    std::vector<size_t> highs = {3, 5, 8, 10};

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
        bool flag = check_result("BC Sort Corner Test", test, res);
        if(!flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }
    

    return 0;
}

int bc_sort_multiple_times(oc::CLP& cmd){

    BASIC_TEST_INIT

    // this test can not pass
    if(role == 0){
        debug_info("RUN BC Sort Multiple Times TEST");
    }

    size_t test_size = 50;
    size_t half_size = test_size / 2;

    // prepare the data.
    aby3::i64Matrix data(test_size, 1);
    aby3::i64Matrix res(test_size, 1);
    for(size_t i=0; i<test_size; i++) {
        data(i, 0) = test_size - 1 - i;
        // res(i, 0) = i;
        if(i < half_size){
            res(i, 0) = half_size + i;
        }
        else{
            res(i, 0) = i - half_size;
        }
    }

    // sort the data.
    std::vector<size_t> lows = {0, half_size};
    std::vector<size_t> highs = {half_size, test_size};

    aby3::sbMatrix enc_data(test_size, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(test_size);
    for(size_t i=0; i<test_size; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    bc_sort_different(enc_data_vec, lows, highs, role, enc, eval, runtime, 625);

    aby3::sbMatrix enc_test(test_size, 64);
    for(size_t i=0; i<test_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    // check the result.
    aby3::i64Matrix test(test_size, 1);
    enc.revealAll(runtime, enc_test, test).get();
    if(role == 0){
        bool flag = check_result("BC Sort Multiple Times Test", test, res);
        if(!flag){
            debug_info("test res = ");
            debug_output_matrix(test);
            debug_info("reference res = ");
            debug_output_matrix(res);
        }
    }

    return 0;
}

int quick_sort_test(oc::CLP &cmd){
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
        debug_info("RUN Quick Sort TEST");
    }

    // setup communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    basic_setup((u64)role, ios, enc, eval, runtime);

    // prepare the data.
    size_t data_size = 1 << 15;
    aby3::i64Matrix data_plain(data_size, 1);
    aby3::i64Matrix data_res(data_size, 1);

    for(size_t i=0; i<data_size; i++){
        data_plain(i, 0) = data_size - i;
        data_res(i, 0) = i + 1;
    }

    aby3::sbMatrix enc_data(data_size, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data_plain, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(data_size);
    for(size_t i=0; i<data_size; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }

    // sort the data.
    size_t min_size = (1 << 5);
    quick_sort_different(enc_data_vec, role, enc, eval, runtime, min_size);

    aby3::sbMatrix enc_test(data_size, 64);
    for(size_t i=0; i<data_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    aby3::i64Matrix test(data_size, 1);
    enc.revealAll(runtime, enc_test, test).get();

    if(role == 0){
        check_result("Quick Sort Test", test, data_res);
    }

    return 0;
}

int quick_sort_with_duplicate_elements_test(oc::CLP& cmd){

    BASIC_TEST_INIT

    if(role == 0){
        debug_info("RUN Quick Sort With Duplicate Elements TEST");
    }

    // prepare the data.
    size_t data_size = 1 << 20;
    size_t bin_size = 1 << 5;
    aby3::i64Matrix data_plain(data_size, 1);
    aby3::i64Matrix data_res(data_size, 1);
    std::vector<int> data_res_vec(data_size);

    for(size_t i=0; i<data_size; i++){
        data_plain(i, 0) = (data_size - i) % bin_size;
        data_res_vec[i] = (i+1) % bin_size;
    }

    std::sort(data_res_vec.begin(), data_res_vec.end());
    for(size_t i=0; i<data_size; i++){
        data_res(i, 0) = data_res_vec[i];
    }

    // todo - add the deplication tags.
    aby3::sbMatrix enc_data(data_size, 64);
    if(role == 0){
        enc.localBinMatrix(runtime, data_plain, enc_data).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_data).get();
    }

    std::vector<aby3::sbMatrix> enc_data_vec(data_size);
    for(size_t i=0; i<data_size; i++){
        aby3::sbMatrix tmp(1, 64);
        tmp.mShares[0](0, 0) = enc_data.mShares[0](i, 0);
        tmp.mShares[1](0, 0) = enc_data.mShares[1](i, 0);
        enc_data_vec[i] = tmp;
    }
    size_t min_size = (1 << 3);

    // tag append, assume that max(data) + log(n) << 64bits, otherwise need to change the tag size.
    quick_sort(enc_data_vec, role, enc, eval, runtime, min_size);
    // tag_append(role, enc_data_vec);


    aby3::sbMatrix enc_test(data_size, 64);
    for(size_t i=0; i<data_size; i++){
        enc_test.mShares[0](i, 0) = enc_data_vec[i].mShares[0](0, 0);
        enc_test.mShares[1](i, 0) = enc_data_vec[i].mShares[1](0, 0);
    }

    aby3::i64Matrix test(data_size, 1);
    enc.revealAll(runtime, enc_test, test).get();

    // if(role == 0){
    //     debug_info("sorted result: ");
    //     debug_output_matrix(test);
    // }

    // if(role == 0){
    //     debug_info("target result: ");
    //     debug_output_matrix(data_res);
    // }

    if(role == 0){
        check_result("Quick Sort with Duplicated Elements Test", test, data_res);
    }
    return 0;
}