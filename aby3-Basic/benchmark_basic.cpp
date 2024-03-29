#include "benchmark_basic.h"
#include "../aby3-RTR/GASTest.h"
#include <cmath>

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
    if(role == 0){
        enc.localBinMatrix(runtime, input_x, enc_x_mat).get();
    }else{
        enc.remoteBinMatrix(runtime, enc_x_mat).get();
    }
    for(size_t i=0; i<data_len; i++){
        enc_x[i].resize(1, 64);
        enc_x[i].mShares[0](0, 0) = enc_x_mat.mShares[0](i, 0);
        enc_x[i].mShares[1](0, 0) = enc_x_mat.mShares[1](i, 0);
    }
    timer.end("data_load");

    // evaluate the sqrt-oram.
    timer.start("oram_construction");
    ABY3SqrtOram oram(data_len, stash_size, pack_size, role, enc, eval, runtime);
    oram.initiate(enc_x);
    timer.end("oram_construction");

    // evaluate the access time.
    std::vector<sbMatrix> access_res(n);
    timer.start("oram_access");
    for(size_t i=0; i<n; i++){
        boolIndex logical_index = boolIndex(i, role);
        access_res[i] = oram.access(logical_index);
    }
    timer.end("oram_access");

    if(role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        stream.close();
    }
    
    return 0;
}