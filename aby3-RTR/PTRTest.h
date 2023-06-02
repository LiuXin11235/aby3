#pragma once
#include "PTRFunction.h"
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>


int test_cipher_index_ptr(oc::CLP& cmd, int n, int m);

int test_cipher_index_ptr_mpi(oc::CLP& cmd, int n, int m, int task_num, int opt_B);

int test_cipher_search_ptr_mpi(oc::CLP& cmd, int n, int m, int task_num, int opt_B);

int test_cipher_rank_ptr_mpi(oc::CLP& cmd, int n, int task_num, int opt_B);

int profile_index(oc::CLP& cmd, int n, int m, int vector_size, int task_num);

int test_vectorization(oc::CLP& cmd, int n, int task_num);