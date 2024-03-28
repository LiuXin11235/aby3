#include "PtAProfile.h"
#include "aby3-Basic/timer.h"
#include "GASTest.h"

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