#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

#include "../aby3-RTR/debug.h"

// #define SHOW_TEST_CASE

int arith_basic_test(oc::CLP& cmd);
int bool_basic_test(oc::CLP& cmd);
int bool_basic_test2(oc::CLP& cmd);
int bool_aggregation_test(oc::CLP& cmd);
int get_first_zero_test(oc::CLP& cmd);

int initialization_test(oc::CLP& cmd);
int shuffle_test(oc::CLP& cmd);
int correlation_test(oc::CLP& cmd);
int communication_test(oc::CLP& cmd);

int pos_map_test(oc::CLP& cmd);
int sqrt_oram_test(oc::CLP& cmd);

int graph_loading_test(oc::CLP& cmd);
int graph_block_fetch_test(oc::CLP& cmd);
int basic_graph_query_test(oc::CLP& cmd);

bool check_result(const std::string& func_name, aby3::i64Matrix& test,
                  aby3::i64Matrix& res);

bool check_result(const std::string& func_name, std::vector<aby3::i64Matrix> test, std::vector<aby3::i64Matrix> res);

template <aby3::Decimal D>
bool check_result(const std::string& func_name, aby3::f64Matrix<D>& test,
                  aby3::f64Matrix<D>& res) {
    return check_result(func_name, test.i64Cast(), res.i64Cast());
}

bool check_result(const std::string& func_name, aby3::i64 test, aby3::i64 res);

template <typename T>
typename std::enable_if<std::is_pod<T>::value, bool>::type check_result(
    const std::string& func_name, std::vector<T>& test, std::vector<T>& res) {
    bool check_flag = true;
    for (size_t i = 0; i < test.size(); i++) {
        if (test[i] != res[i]) {
            check_flag = false;
        }
    }
    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
#ifdef SHOW_TEST_CASE
        debug_info("test case: ");
        debug_output_vector(test);
        debug_info("expected result: ");
        debug_output_vector(res);
#endif
    } else {
        debug_info("\033[32m" + func_name + " SUCCESS!" + "\033[0m\n");
    }
    return check_flag;
}

template <typename T>
typename std::enable_if<std::is_pod<T>::value, bool>::type check_result(
    const std::string& func_name, T& test, T& res) {
    bool check_flag = (test == res);
    if (!check_flag) {
        debug_info("\033[31m" + func_name + " ERROR !" + "\033[0m\n");
#ifdef SHOW_TEST_CASE
        debug_info("test case: ");
        debug_output_value(test);
        debug_info("expected result: ");
        debug_output_value(res);
#endif
    }
    return check_flag;
}