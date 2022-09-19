#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>


int test_mul(oc::CLP& cmd);

int test_gt(oc::CLP& cmd);

// int test_eq(oc::CLP& cmd);

int test_argsort(oc::CLP& cmd, int rtrFlag);

int basic_performance(oc::CLP& cmd, int n, int repeats, std::map<std::string, std::vector<double>>& dict);