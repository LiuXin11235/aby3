#include "Test.h"

#include <chrono>
#include <random>
#include <thread>
#include <math.h>
#include <vector>
#include <algorithm>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/CLP.h>

#include "../aby3-GORAM-Core/Basics.h"
#include "../aby3-GORAM-Core/Sort.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-ML/main-logistic.h"
#include "../aby3-ML/Regression.h"  

using namespace oc;
using namespace aby3;
using namespace Eigen;

#define P0_IP "10.233.105.209" //"127.0.0.1"
#define P1_IP "10.233.105.219" //"127.0.0.1"

#define PORT10 1412
#define PORT20 1210
#define PORT21 1211

int lr_test(oc::CLP& cmd){

    int role = -1;
    if (cmd.isSet("role")) {
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if (role == -1) {
        throw std::runtime_error(LOCATION);
    }
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;

    if(role == 0){
        debug_info("RUN Logistic Regression TEST");
    }

    // prepare the data.
    // int N = 1000, D = 100, B = 128, IT = 1000;
    // int testN = 1000;
    int B = 128, IT = 100; // batch size and number of iterations!!
    int N = 100, testN = 100, D = 100;

    RegressionParam params;
    params.mBatchSize = B;
    params.mIterations = IT;
    params.mLearningRate = 1.0 / (1 << 3);

    auto next = (role + 1) % 3;
    auto prev = (role + 2) % 3;
    auto cNameNext = std::to_string(std::min(role, next)) + std::to_string(std::max(role, next));
    auto cNamePrev = std::to_string(std::min(role, prev)) + std::to_string(std::max(role, prev));

    auto modeNext = role < next ? SessionMode::Server : SessionMode::Client;
    auto modePrev = role < prev ? SessionMode::Server : SessionMode::Client;

    // // The single-machine version: original
    // auto portNext = 1212 + std::min(role, next);
    // auto portPrev = 1212 + std::min(role, prev);
    // Session epNext(ios, "127.0.0.1", portNext, modeNext, cNameNext);
    // Session epPrev(ios, "127.0.0.1", portPrev, modePrev, cNamePrev);    
    
    // Initialize sessions based on role in distributed version
    std::string IPNext;
    std::string IPPrev;    
    std::string portNext;
    std::string portPrev;
    if (role == 0) {
        IPNext = std::string(P0_IP);
        IPPrev = std::string(P0_IP);
        portNext = std::to_string(PORT10); // 0 as server to receive , 1 as client to send
        portPrev = std::to_string(PORT20); // 0 as server, 2 as client
    } else if (role == 1) {
        IPNext = std::string(P1_IP);
        IPPrev = std::string(P0_IP);
        portNext = std::to_string(PORT21); // 1 as server, 2 as client
        portPrev = std::to_string(PORT10); // 0 as server, 1 as client
    } else {
        IPNext = std::string(P0_IP);
        IPPrev = std::string(P1_IP);
        portNext = std::to_string(PORT20); // 0 as server, 2 as client
        portPrev = std::to_string(PORT21); // 1 as server, 2 as client
    }
    Session epNext(ios, IPNext + ":" + portNext, modeNext, cNameNext);
    Session epPrev(ios, IPPrev + ":" + portPrev, modePrev, cNamePrev);

    // Initialize channels
    auto chlNext = epNext.addChannel();
    auto chlPrev = epPrev.addChannel();

    chlNext.waitForConnection();
    chlPrev.waitForConnection();

    // chlNext.send(role);
    // chlPrev.send(role);
    // u64 prevAct, nextAct;
    // chlNext.recv(nextAct);
    // chlPrev.recv(prevAct);

    if(cmd.isSet("LR_train")){
      aby3::logistic_main_3pc_sh(N, D, B, IT, testN, role, true, cmd, epPrev, epNext);
    }
    if(cmd.isSet("LR_test")){
      aby3::logistic_main_3pc_sh_test(N, D, B, IT, testN, role, true, cmd, epPrev, epNext);
    }
    // logistic_main_3pc_sh(cmd);

    return 0;
}
