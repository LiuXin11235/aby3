#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include "PtATasks.h"

#define SETUP_PROCESS \
    int rank, size;  \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  \
    MPI_Comm_size(MPI_COMM_WORLD, &size); \
    int role = -1; \
    if (cmd.isSet("role")) { \
        auto keys = cmd.getMany<int>("role"); \
        role = keys[0]; \
    } \
    if (role == -1) { \
        throw std::runtime_error(LOCATION); \
    } \
    IOService ios; \
    Sh3Encryptor enc; \
    Sh3Evaluator eval; \
    Sh3Runtime runtime; \
    multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
    

int test_cipher_index_pta(oc::CLP& cmd);