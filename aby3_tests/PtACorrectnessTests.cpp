#include "Test.h"
#include <mpi.h>

#include "../aby3-RTR/PtATests.h"
#include "../aby3-RTR/GASTest.h"
#include "../aby3-Basic/timer.h"
#include "../aby3-Basic/Basics.h"

using namespace oc;
using namespace aby3;


int correctness_cipher_index_pta(oc::CLP& cmd){
    size_t n = 1, m=1<<10, optB=128;
    int task_num;
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);

    SETUP_PROCESS
    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, int, si64, si64, CipherIndex>(task_num, optB, role, enc, runtime, eval);
    ptaTask->circuit_construct({n}, {m});
    ptaTask->set_default_value(GET_ZERO_SHARE);

    // data loading.
    std::vector<si64> inputX(n); std::vector<int> inputY; 
    std::vector<si64> inputV; std::vector<si64> res(n);

    // generate the test data.
    size_t partial_len = ptaTask->subTask->get_partial_m_lens();
    i64Matrix pinput_x(n, 1); i64Matrix pinput_v(partial_len, 1);
    pinput_x(0, 0) = 1;
    inputY.resize(partial_len); inputV.resize(partial_len);
    for(size_t i=0; i<partial_len; i++){
        inputY[i] = i + ptaTask->m_start;
        pinput_v(i, 0) = i + ptaTask->m_start;
    }

    aby3::si64Matrix sinputX(n, 1);
    aby3::si64Matrix sinputV(partial_len, 1);
    if(role == 0){
        enc.localIntMatrix(runtime, pinput_x, sinputX).get();
        enc.localIntMatrix(runtime, pinput_v, sinputV).get();
    }
    else{
        enc.remoteIntMatrix(runtime, sinputX).get();
        enc.remoteIntMatrix(runtime, sinputV).get();
    }
    inputX[0].mData[0] = sinputX.mShares[0](0, 0); 
    inputX[0].mData[1] = sinputX.mShares[1](0, 0);
    
    res[0] = GET_ZERO_SHARE;
    for(size_t i=0; i<partial_len; i++){
        // inputV[i] = sinputV(i, 0);
        inputV[i].mData[0] = sinputV.mShares[0](i, 0);
        inputV[i].mData[1] = sinputV.mShares[1](i, 0);
    }

    std::string str_inputx = "x = ";
    std::string str_inputv = "v = ";
    
    aby3::i64Matrix test_plain_x = back2plain(role, inputX, enc, eval, runtime);
    aby3::i64Matrix test_plain_v = back2plain(role, inputV, enc, eval, runtime);
    str_inputx += std::to_string(test_plain_x(0, 0));
    for(size_t i=0; i<partial_len; i++){
        str_inputv += std::to_string(test_plain_v(i, 0)) + " ";
    }

    debug_mpi(rank, role, str_inputx);
    debug_mpi(rank, role, str_inputv);

    // evaluate the task.
    ptaTask->set_selective_value(inputV.data(), 0);
    ptaTask->circuit_evaluate(inputX.data(), inputY.data(), inputV.data(), res.data());

    // check the result.   
    aby3::i64Matrix test_res(1, 1);
    aby3::si64Matrix runtime_res(1, 1);
    // runtime_res(0, 0) = res[0];
    runtime_res.mShares[0](0, 0) = res[0].mData[0];
    runtime_res.mShares[1](0, 0) = res[0].mData[1];
    enc.revealAll(runtime, runtime_res, test_res).get();

    if(rank == 0 && role == 0){
        check_result("pta cipher_index", test_res, pinput_x);
    }
    return 0;
}