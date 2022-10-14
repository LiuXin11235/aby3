#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

int dis_test_mul(oc::CLP& cmd);

// performance test.
int dis_basic_performance(oc::CLP& cmd, int n, int repeats, std::map<std::string, double>& dict);