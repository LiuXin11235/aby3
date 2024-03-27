#include "PtATests.h"

int test_cipher_index_pta(oc::CLP& cmd, size_t n, size_t m, int task_num, int opt_B){

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
    multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
    
}