#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include "../aby3-RTR/debug.h"


int arith_basic_test(oc::CLP& cmd);
int bool_basic_test(oc::CLP& cmd);
int initialization_test(oc::CLP& cmd);
int shuffle_test(oc::CLP& cmd);
int correlation_test(oc::CLP& cmd);

bool check_result(const std::string& func_name, aby3::i64Matrix& test, aby3::i64Matrix& res);

template <aby3::Decimal D>
bool check_result(const std::string& func_name, aby3::f64Matrix<D>& test, aby3::f64Matrix<D>& res){
    return check_result(func_name, test.i64Cast(), res.i64Cast());
}

bool check_result(const std::string& func_name, aby3::i64 test, aby3::i64 res);

template<typename T>
typename std::enable_if<std::is_pod<T>::value, bool>::type 
check_result(const std::string& func_name, std::vector<T> &test, std::vector<T> &res){
    bool check_flag = true;
    for (size_t i = 0; i < test.size(); i++){
        if (test[i] != res[i]){
            check_flag = false;
        }
    }
    if(!check_flag){
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
    }
    else{
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }
    return check_flag;
}