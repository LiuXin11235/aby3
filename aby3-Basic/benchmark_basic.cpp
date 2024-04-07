#include "benchmark_basic.h"
#include "../aby3-RTR/GASTest.h"
#include <cmath>

static size_t MAX_ENC_SIZE = 1 << 25;

using namespace oc;
using namespace aby3;

int sqrt_oram_benchmark(oc::CLP& cmd){
    size_t m, n; 
    get_value("N", cmd, n);
    get_value("M", cmd, m);
    size_t data_len = m;
    std::string logging_file;
    get_value("logFile", cmd, logging_file);

    Timer& timer = Timer::getInstance();
    timer.start("time_setup");
    SETUP_SINGLE_PROCESS
    timer.end("time_setup");

    timer.start("data_load");
    // generate_data, stash_size = sqrt(n)
    size_t stash_size = (size_t) sqrt((double)data_len);
    size_t pack_size = 8;
    size_t block_size = 1;

    // construct the sqrt-oram.
    // std::vector<i64Matrix> input_x(n);
    i64Matrix input_x;
    sbMatrix enc_x_mat;
    std::vector<sbMatrix> enc_x(data_len);

    input_x.resize(data_len, 1);
    for(size_t i=0; i<data_len; i++){
        input_x(i, 0) = i;
    }
    enc_x_mat.resize(data_len, 64);

    size_t round = (size_t) ceil(data_len / (double)MAX_ENC_SIZE);
    size_t last_len = data_len - (round - 1) * MAX_ENC_SIZE;
    for(size_t i=0; i<round; i++){
        size_t enc_len = (i == round - 1) ? last_len : MAX_ENC_SIZE;
        i64Matrix x_block = input_x.block(i * MAX_ENC_SIZE, 0, enc_len, 1);
        sbMatrix _enc_block(enc_len, 64);
        if(role == 0){
            enc.localBinMatrix(runtime, x_block, _enc_block).get();
        }else{
            enc.remoteBinMatrix(runtime, _enc_block).get();
        }
        // enc_x_mat.block(i * MAX_ENC_SIZE, 0, enc_len, 64) = _enc_block;
        for(size_t j=0; j<enc_len; j++){
            enc_x_mat.mShares[0](i * MAX_ENC_SIZE + j, 0) = _enc_block.mShares[0](j, 0);
            enc_x_mat.mShares[1](i * MAX_ENC_SIZE + j, 0) = _enc_block.mShares[1](j, 0);
        }
    }

    for(size_t i=0; i<data_len; i++){
        enc_x[i].resize(1, 64);
        enc_x[i].mShares[0](0, 0) = enc_x_mat.mShares[0](i, 0);
        enc_x[i].mShares[1](0, 0) = enc_x_mat.mShares[1](i, 0);
    }
    timer.end("data_load");

    std::cout << "after data loading" << std::endl;

    // evaluate the sqrt-oram.
    timer.start("oram_construction");
    ABY3SqrtOram oram(data_len, stash_size, pack_size, role, enc, eval, runtime);
    oram.initiate(enc_x);
    timer.end("oram_construction");

    std::cout << "after oram construction" << std::endl;

    // evaluate the access time.
    std::vector<sbMatrix> access_res(n);
    timer.start("oram_access");
    for(size_t i=0; i<n; i++){
        boolIndex logical_index = boolIndex(i, role);
        access_res[i] = oram.access(logical_index);
    }
    timer.end("oram_access");

    std::cout << "after oram access" << std::endl;

    if(role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        stream.close();
    }
    
    return 0;
}