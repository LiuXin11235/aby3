#include "PtATests.h"
#include "GASTest.h"
#include "../aby3-Basic/timer.h"

using namespace oc;
using namespace aby3;


int test_cipher_index_pta(oc::CLP& cmd){

    size_t n, m, optB;
    int task_num;
    std::string logging_file;
    get_value("N", cmd, n); get_value("M", cmd, m);
    get_value("B", cmd, optB);
    get_value("logFile", cmd, logging_file);
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);

    Timer& timer = Timer::getInstance();

    timer.start("time_setup");
    SETUP_PROCESS
    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, int, si64, si64, CipherIndex>(task_num, optB, role, enc, runtime, eval);
    ptaTask->circuit_construct({n}, {m});
    ptaTask->set_default_value(GET_ZERO_SHARE);
    timer.end("time_setup");

    timer.start("time_data_prepare");
    // data loading.
    std::vector<si64> inputX; std::vector<int> inputY; 
    std::vector<si64> inputV; std::vector<si64> res(n);
    std::tie(inputX, inputY, inputV) = ptaTask->subTask->data_loading();
    // pass the seletive value.
    ptaTask->set_selective_value(inputV.data(), 0);
    timer.end("time_data_prepare");
    
    // evaluate the task.
    timer.start("time_circuit_evaluate");
    ptaTask->circuit_evaluate(inputX.data(), inputY.data(), inputV.data(), res.data());
    timer.end("time_circuit_evaluate");

    // print the time.
    if(rank == 0 && role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        ptaTask->print_time_profile(stream);
        stream.close();
    }

    return 0;
}

int test_max_pta(oc::CLP& cmd){

    size_t n, m, optB;
    int task_num;
    std::string logging_file;
    get_value("N", cmd, n); get_value("M", cmd, m);
    get_value("B", cmd, optB);
    get_value("logFile", cmd, logging_file);
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);

    if(n != 1) n = 1;

    Timer& timer = Timer::getInstance();

    timer.start("time_setup");
    SETUP_PROCESS
    // construct the task.
    auto ptaTask = new ABY3MPITask<int, si64, si64, si64, Max>(task_num, optB, role, enc, runtime, eval);
    ptaTask->circuit_construct({n}, {m});
    ptaTask->set_default_value(GET_ZERO_SHARE);
    timer.end("time_setup");

    timer.start("time_data_prepare");
    // data loading.
    std::vector<si64> inputY = ptaTask->subTask->data_loading();
    std::vector<si64> res(n);
    std::vector<int> inputX(1);
    timer.end("time_data_prepare");

    // evaluate the task.
    timer.start("time_circuit_evaluate");
    ptaTask->circuit_evaluate(inputX.data(), inputY.data(), nullptr, res.data());
    timer.end("time_circuit_evaluate");

    // print the time.
    if(rank == 0 && role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        ptaTask->print_time_profile(stream);
        stream.close();
    }

    return 0;
}