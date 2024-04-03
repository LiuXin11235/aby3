#include "PtATests.h"
#include "GASTest.h"
#include "../aby3-Basic/timer.h"
#include "debug.h"
// #include "../aby3_tests/Test.h"

using namespace oc;
using namespace aby3;

aby3::si64 get_share(int pIdx, int target_val){
    aby3::si64 dval;
    switch(pIdx){
        case 0:
            dval.mData[0] = 0;
            dval.mData[1] = 0;
            break;
        case 1:
            dval.mData[0] = target_val;
            dval.mData[1] = 0;
            break;
        case 2:
            dval.mData[0] = 0;
            dval.mData[1] = target_val;
            break;
        default:
            THROW_RUNTIME_ERROR("Invalid party index - " + std::to_string(pIdx));
    }
    return dval;
}

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
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->circuit_construct({n}, {m});
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
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->circuit_construct({n}, {m});
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

int test_sort_pta(oc::CLP& cmd){
    size_t n, m, optB;
    int task_num;
    std::string logging_file;
    get_value("N", cmd, n); get_value("M", cmd, m);
    get_value("B", cmd, optB);
    get_value("logFile", cmd, logging_file);
    MPI_Comm_size(MPI_COMM_WORLD, &task_num);

    if(n != m){
        THROW_RUNTIME_ERROR("n should be equal to m for sort test!");
    };

    Timer& timer = Timer::getInstance();

    timer.start("time_setup");
    SETUP_PROCESS
    // construct the task.
    auto ptaTask_step1 = new ABY3MPITask<si64, si64, si64, si64, Rank>(task_num, optB, role, enc, runtime, eval);
    ptaTask_step1->set_default_value(GET_ZERO_SHARE);
    ptaTask_step1->circuit_construct({n}, {m});
    auto ptaTask_step2 = new ABY3MPITask<int, si64, si64, si64, CipherIndex>(task_num, optB, role, enc, runtime, eval);
    ptaTask_step2->set_default_value(GET_ZERO_SHARE);
    ptaTask_step2->circuit_construct({n}, {m});
    timer.end("time_setup");

    timer.start("time_data_prepare");
    // data loading.
    std::vector<si64> inputY; std::vector<si64> inputX;
    std::tie(inputX, inputY) = ptaTask_step1->subTask->data_loading();
    std::vector<si64> res(n);
    std::vector<si64> sort_res(n);
    ptaTask_step2->set_selective_value(inputY.data(), 0);
    std::vector<int> range_index(n);
    for(size_t i=0; i<n; i++) range_index[i] = i;
    timer.end("time_data_prepare");

    // evaluate the task (1 step)
    timer.start("time_circuit_evaluate_step1");
    ptaTask_step1->circuit_evaluate(inputX.data(), inputY.data(), nullptr, res.data());
    timer.end("time_circuit_evaluate_step1");
    MPI_Barrier(MPI_COMM_WORLD);
    
    // resharing the result and prepare for the next step.
    timer.start("time_data_prepare_step2");
    size_t sharing_length;
    aby3::si64* partial_data;
    if(rank == 0){
        sharing_length = ptaTask_step1->m;
        partial_data = res.data();
    }
    else{
        sharing_length = (ptaTask_step1->m_end - ptaTask_step1->m_start) + 1;
        partial_data = new aby3::si64[sharing_length];
    }
    ptaTask_step2->data_sharing<aby3::si64>(partial_data, sharing_length, 1);
    timer.end("time_data_prepare_step2");
    MPI_Barrier(MPI_COMM_WORLD);

    timer.start("time_circuit_evaluate_step2");
    ptaTask_step2->circuit_evaluate(range_index.data(), partial_data, inputY.data(), sort_res.data());
    timer.end("time_circuit_evaluate_step2");

    if(rank == 0 && role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        timer.print_total("milliseconds", stream);
        stream << "step1: " << std::endl;
        ptaTask_step1->print_time_profile(stream);
        stream << "step2: " << std::endl;
        ptaTask_step2->print_time_profile(stream);
    }

    return 0;
}

int test_sum_pta(oc::CLP& cmd){

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
    auto ptaTask = new ABY3MPITask<int, si64, si64, si64, Sum>(task_num, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);
    ptaTask->circuit_construct({n}, {m});
    timer.end("time_setup");

    timer.start("time_data_prepare");
    // data loading.
    std::vector<si64> partialX = ptaTask->subTask->data_loading();
    std::vector<si64> res(n);
    std::vector<int> inputX(1);
    timer.end("time_data_prepare");

    // evaluate the task.
    timer.start("time_circuit_evaluate");
    ptaTask->circuit_evaluate(inputX.data(), partialX.data(), nullptr, res.data());
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

int test_metric_pta(oc::CLP& cmd){

    size_t n, m, optB;
    size_t k = 1;
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
    auto ptaTask = new ABY3MPITask<std::vector<si64>, std::vector<si64>, si64, si64, BioMetric>(task_num, optB, role, enc, runtime, eval);
    ptaTask->set_default_value(get_share(role, m+1));
    ptaTask->circuit_construct({n}, {m});
    timer.end("time_setup");

    timer.start("time_data_prepare");
    // data loading.
    std::vector<std::vector<si64>> inputX; std::vector<std::vector<si64>> inputY;
    std::tie(inputX, inputY) = ptaTask->subTask->data_loading();
    std::vector<si64> res(n);
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