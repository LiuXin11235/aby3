#include "PtAProfile.h"
#include "aby3-Basic/timer.h"
#include "GASTest.h"
#include "BuildingBlocks.h"

using namespace oc;
using namespace aby3;

// #define SEE_VECTOR

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

int communication_profile(oc::CLP& cmd){

    Timer& timer = Timer::getInstance();

    timer.start("time_setup");
    SETUP_PROCESS
    timer.end("time_setup");

    // test the communication & bandwidth.
    std::string logging_file;
    get_value("logFile", cmd, logging_file);
    size_t start_b, end_b;
    get_value("startB", cmd, start_b);
    get_value("endB", cmd, end_b);

    size_t b = start_b;

    aby3::si64Matrix small_data(1, 1);
    init_ones(role, enc, runtime, small_data, 1);
    MPI_Barrier(MPI_COMM_WORLD);

    timer.start("time_latency");
    runtime.mComm.mNext.asyncSend(small_data.mShares[0].data(), small_data.mShares[0].size());
    auto fu = runtime.mComm.mPrev.asyncRecv(small_data.mShares[1].data(), small_data.mShares[1].size());
    fu.get();
    MPI_Barrier(MPI_COMM_WORLD);
    timer.end("time_latency");


    while(b <= end_b){

        aby3::si64Matrix data(b, 1);
        init_zeros(role, enc, runtime, data, b);

        MPI_Barrier(MPI_COMM_WORLD);

        std::string comm_key = "time_communication." + std::to_string(b);
        timer.start(comm_key);

        auto fu_send = runtime.mComm.mNext.asyncSendFuture(data.mShares[0].data(), data.mShares[0].size());
        auto fu_recv = runtime.mComm.mPrev.asyncRecv(data.mShares[1].data(), data.mShares[1].size());
        fu_send.get();
        fu_recv.get();

        timer.end(comm_key);

        b *= 4;
    }

    if(rank ==0 && role ==0){
        std::ofstream stream(logging_file, std::ios::app);
        stream << "tasknum: " << std::to_string(size) << std::endl;
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
    else if(target_task == "sort"){
        sort_profile(cmd);
        return 0;
    }
    else if(target_task == "sum"){
        sum_profile(cmd);
        return 0;
    }
    else if(target_task == "metric"){
        metric_profile(cmd);
        return 0;
    }
    else{
        throw std::runtime_error(LOCATION);
    }
    return 0;
}


std::pair<size_t, double> get_optimal_vector_size(size_t b_start, size_t b_end, size_t gap, std::function<std::pair<double, double>(size_t)> evaluate_task, std::string logging_file){
    size_t n=1; 
    double last_ratio = -1;
    double ratio = -1;
    size_t b = b_start;
    double time_c = -1;
    double time_c_all = 0;
    size_t repeat_times = 5;

    std::ofstream logging_stream;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        logging_stream.open(logging_file, std::ios::app);
        logging_stream << "start_b: " << b_start << std::endl;
        logging_stream << "end_b: " << b_end << std::endl;
        logging_stream << "gap: " << gap << std::endl;

        logging_stream << "exponential probing" << std::endl;
    }

    // exponentially increase the batch size for probing.
    while(b <= b_end){

        time_c_all = 0;
        for(size_t i=0; i<repeat_times; i++){
            std::tie(time_c, ratio) = evaluate_task(b);
            time_c_all += time_c;
        }
        time_c = time_c_all / repeat_times;
        ratio = time_c / b;

#ifndef SEE_VECTOR
        if(last_ratio > 0){
            if(ratio > last_ratio * 1.1){
                if(rank == 0){
                    logging_stream << "b: " << b << std::endl;
                    logging_stream << "time_c: " << time_c << std::endl;
                    logging_stream << "ratio: " << ratio << std::endl;
                }
                break;
            }
        }
#endif
        last_ratio = ratio;

        if(rank == 0){
            logging_stream << "b: " << b << std::endl;
            logging_stream << "time_c: " << time_c << std::endl;
            logging_stream << "ratio: " << ratio << std::endl;
        }

        b *= 2;
    }

#ifdef SEE_VECTOR
    return {1, 0};
#endif

    return {b/2, time_c};

    // then binary search the optimal batch size.
    size_t binary_start_b = b / 2;
    size_t binary_end_b = b;

    if(rank == 0){
        logging_stream << "binary probing" << std::endl;
    }

    while(binary_end_b - binary_start_b > gap){
        b = (binary_start_b + binary_end_b) / 2;
        time_c_all = 0;

        for(size_t i=0; i<repeat_times; i++){
            std::tie(time_c, ratio) = evaluate_task(b);
            time_c_all += time_c;
        }

        time_c = time_c_all / repeat_times;
        ratio = time_c / b;

        if(ratio > last_ratio){
            binary_end_b = b;
            if(rank == 0){
                logging_stream << "b: " << b << std::endl;
                logging_stream << "time_c: " << time_c << std::endl;
                logging_stream << "ratio: " << ratio << std::endl;
            }
        }
        else break;
    }

    // get the final result.
    size_t optimal_vector_size = (binary_start_b + binary_end_b) / 2;
    std::tie(time_c, ratio) = evaluate_task(optimal_vector_size);

    return {optimal_vector_size, time_c};
}

double get_unit_time(size_t b, std::function<std::pair<double, double>(size_t)> evaluate_task, std::string logging_file){

    double time_c = -1;
    double ratio = -1;
    double time_c_all = 0;
    size_t repeat_times = 5;

    for(size_t i=0; i<repeat_times; i++){
        std::tie(time_c, ratio) = evaluate_task(b);
        time_c_all += time_c;
    }
    time_c = time_c_all / repeat_times;

    return time_c;
}

int cipher_index_profile(oc::CLP& cmd){

    // setup the process.
    SETUP_PROCESS

    // prepare the profiler.
    PROFILER_PREPARE

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

        double _time_c = timer.get_time("time_circuit_evaluate", "milliseconds");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();
        MPI_Barrier(MPI_COMM_WORLD);

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}

int max_profile(oc::CLP& cmd){

    // setup the process.
    SETUP_PROCESS

    // prepare the profiler.
    PROFILER_PREPARE

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

        double _time_c = timer.get_time("time_circuit_evaluate", "milliseconds");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();
        MPI_Barrier(MPI_COMM_WORLD);

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}


int sort_profile(oc::CLP& cmd){

    // setup the process.
    SETUP_PROCESS    

    // prepare the profiler.
    PROFILER_PREPARE

    // construct the task.
    auto ptaTask_step1 = new ABY3MPITask<si64, si64, si64, si64, Rank>(size, start_b, role, enc, runtime, eval);
    auto ptaTask_step2 = new ABY3MPITask<int, si64, si64, si64, CipherIndex>(size, start_b, role, enc, runtime, eval);
    ptaTask_step1->set_default_value(GET_ZERO_SHARE);
    ptaTask_step2->set_default_value(GET_ZERO_SHARE);

    // define the task evaluation function.
    auto evaluate_task = [&](size_t b){

        // get the timer.
        Timer& timer = Timer::getInstance();
        size_t m= (size_t) sqrt(size * b);

        ptaTask_step1->optimal_block = b;
        ptaTask_step1->circuit_construct({m}, {m});
        ptaTask_step2->optimal_block = b;
        ptaTask_step2->circuit_construct({m}, {m});

        // data loading.
        std::vector<si64> inputX; std::vector<si64> inputY;
        std::tie(inputX, inputY) = ptaTask_step1->subTask->data_loading();
        std::vector<int> range_index(n);
        for(size_t i=0; i<n; i++) range_index[i] = i;   

        // pass the seletive value.
        auto fakeX = ptaTask_step1->fake_repeatX(inputX.data());
        auto fakeY = ptaTask_step1->fake_repeatY(inputY.data());
        auto fakeX2 = ptaTask_step2->fake_repeatX(range_index.data());
        ptaTask_step2->set_selective_value(inputY.data(), 0);

        // evaluate the task.
        timer.start("time_circuit_evaluate");
        ptaTask_step1->subTask->circuit_profile(fakeX, fakeY, ptaTask_step1->selectV);
        ptaTask_step2->subTask->circuit_profile(fakeX2, fakeY, ptaTask_step2->selectV);
        timer.end("time_circuit_evaluate");

        double _time_c = timer.get_time("time_circuit_evaluate", "milliseconds");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();
        MPI_Barrier(MPI_COMM_WORLD);

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}

int sum_profile(oc::CLP& cmd){

    // setup the process.
    SETUP_PROCESS

    // prepare the profiler.
    PROFILER_PREPARE

    // construct the task.
    auto ptaTask = new ABY3MPITask<int, si64, si64, si64, Sum>(size, start_b, role, enc, runtime, eval);
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
        std::vector<si64> partialY = ptaTask->subTask->data_loading();
        std::vector<si64> res(n);
        std::vector<int> inputX(1);

        auto fakeY = ptaTask->fake_repeatY(partialY.data());
        auto fakeX = ptaTask->fake_repeatX(inputX.data());

        timer.start("time_circuit_evaluate");
        ptaTask->subTask->circuit_profile(fakeX, fakeY, ptaTask->selectV);
        timer.end("time_circuit_evaluate");

        double _time_c = timer.get_time("time_circuit_evaluate", "milliseconds");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();
        MPI_Barrier(MPI_COMM_WORLD);
        
        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER

    return 0;
}

int metric_profile(oc::CLP& cmd){

    // setup the process.
    SETUP_PROCESS

    // prepare the profiler.
    PROFILER_PREPARE

    // construct the task.
    auto ptaTask = new ABY3MPITask<std::vector<si64>, std::vector<si64>, si64, si64, BioMetric>(size, start_b, role, enc, runtime, eval);
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
        std::vector<std::vector<si64>> inputX; std::vector<std::vector<si64>> inputY;
        std::tie(inputX, inputY) = ptaTask->subTask->data_loading();
        std::vector<si64> res(n);

        auto fakeX = ptaTask->fake_repeatX(inputX.data());
        auto fakeY = ptaTask->fake_repeatY(inputY.data());

        timer.start("time_circuit_evaluate");
        ptaTask->subTask->circuit_profile(fakeX, fakeY, ptaTask->selectV);
        timer.end("time_circuit_evaluate");

        double _time_c = timer.get_time("time_circuit_evaluate", "milliseconds");
        synchronized_time(role, _time_c, runtime);
        MPI_Bcast(&_time_c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double ratio = _time_c / b;
        timer.clear_records();
        MPI_Barrier(MPI_COMM_WORLD);

        return std::make_pair(_time_c, ratio);
    };

    // get the optimal vector size.
    PROFILER_RECORDER
    

    return 0;
}