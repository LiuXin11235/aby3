#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include <cmath>
#include "PTRFunction.h"

using namespace oc;
using namespace aby3;
using namespace std;

int profile_cipher_index(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);

int profile_average(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);

int profile_rank(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);

int profile_sort(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);

int profile_max(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap);

int profile_bio_metric(oc::CLP& cmd, size_t n, size_t m, size_t k, int vector_size_start, double epsilon, size_t gap);

int profile_mean_distance(oc::CLP& cmd, size_t n, size_t m, size_t k, int vector_size_start, double epsilon, size_t gap);