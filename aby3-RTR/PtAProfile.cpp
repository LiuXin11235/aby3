#include "PtAProfile.h"
#include "aby3-Basic/timer.h"
#include "GASTest.h"
#include "BuildingBlocks.h"

using namespace oc;
using namespace aby3;

int pta_system_profile(oc::CLP& cmd){

    Timer& timer = Timer::getInstance();
    timer.start("time_setup");
    SETUP_PROCESS
    timer.end("time_setup");

    std::string logging_file;
    get_value("logFile", cmd, logging_file);

    if(rank == 0 && role == 0){
        std::ofstream stream(logging_file, std::ios::app);
        stream << "task num: " << size << std::endl;
        timer.print_total("milliseconds", stream);
        stream.close();
    }

    return 0;
}


int task_profile(oc::CLP& cmd){
    std::string target_task;
    get_value("task", cmd, target_task);

    if(target_task == "cipher_index"){
        cipher_index_profile(cmd);
        return 0;
    }
    else if(target_task == "max"){
        max_profile(cmd);
        return 0;
    }
    else{
        throw std::runtime_error(LOCATION);
    }
    return 0;
}


std::pair<size_t, double> get_optimal_vector_size(size_t b_start, size_t b_end, size_t gap, std::function<std::pair<double, double>(size_t)> evaluate_task){
    size_t n=1; 
    double last_ratio = -1;
    double ratio = -1;
    size_t b = b_start;
    double time_c = -1;

    // exponentially increase the batch size for probing.
    while(b < b_end){
        if(b > b_end) b = b_end;

        std::tie(ratio, time_c) = evaluate_task(b);

        if(last_ratio > 0){
            if(ratio < last_ratio) break;
        }
        last_ratio = ratio;

        b *= 2;
    }

    // then binary search the optimal batch size.
    size_t binary_start_b = b / 2;
    size_t binary_end_b = b;

    while(binary_end_b - binary_start_b > gap){
        b = (binary_start_b + binary_end_b) / 2;
        std::tie(ratio, time_c) = evaluate_task(b);
        if(ratio > last_ratio) binary_end_b = b;
        else break;
    }

    // get the final result.
    size_t optimal_vector_size = (binary_start_b + binary_end_b) / 2;
    std::tie(ratio, time_c) = evaluate_task(optimal_vector_size);

    return {optimal_vector_size, time_c};
}


int cipher_index_profile(oc::CLP& cmd){

    // prepare the profiler.
    PROFILER_PREPARE

    // setup the process.
    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<si64, int, si64, si64, CipherIndex>(size, start_b, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);

    // define the task evaluation function.
    auto evaluate_task = [&](size_t b){

        // get the timer.
        Timer& timer = Timer::getInstance();
        size_t m=size * b;

        ptaTask->optimal_block = b;
        ptaTask->circuit_construct({n}, {m});

        // data loading.
        std::vector<si64> inputX; std::vector<int> inputY; std::vector<si64> inputV;
        std::tie(inputX, inputY, inputV) = ptaTask->subTask->data_loading();

        // pass the seletive value.
        ptaTask->set_selective_value(inputV.data(), 0);
        auto fakeX = ptaTask->fake_repeatX(inputX.data());
        auto fakeY = ptaTask->fake_repeatY(inputY.data());

        // evaluate the task.
        timer.start("time_circuit_evaluate");
        ptaTask->subTask->circuit_profile(fakeX, fakeY, ptaTask->selectV);
        timer.end("time_circuit_evaluate");

        double _time_c = timer.get_time("time_circuit_evaluate");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}

int max_profile(oc::CLP& cmd){

    // prepare the profiler.
    PROFILER_PREPARE

    // setup the process.
    SETUP_PROCESS

    // construct the task.
    auto ptaTask = new ABY3MPITask<int, si64, si64, si64, Max>(size, start_b, role, enc, runtime, eval);
    ptaTask->set_default_value(GET_ZERO_SHARE);

    // define the task evaluation function.
    auto evaluate_task = [&](size_t b){

        // get the timer.
        Timer& timer = Timer::getInstance();
        size_t m=size * b;

        if(n != 1) n = 1;

        ptaTask->optimal_block = b;
        ptaTask->circuit_construct({n}, {m});

        // data loading.
        std::vector<int> inputX(1);
        std::vector<si64> inputY = ptaTask->subTask->data_loading();

        // pass the seletive value.
        auto fakeX = ptaTask->fake_repeatX(inputX.data());
        auto fakeY = ptaTask->fake_repeatY(inputY.data());

        // evaluate the task.
        timer.start("time_circuit_evaluate");
        ptaTask->subTask->circuit_profile(fakeX, fakeY, ptaTask->selectV);
        timer.end("time_circuit_evaluate");

        double _time_c = timer.get_time("time_circuit_evaluate");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}