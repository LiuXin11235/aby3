#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include "PTRFunction.h"

using namespace oc;
using namespace aby3;
using namespace std;

int profile_cipher_index(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);
