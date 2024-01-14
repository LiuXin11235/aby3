#pragma once
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

int arith_basic_test(oc::CLP& cmd);
int bool_basic_test(oc::CLP& cmd);

bool check_result(const std::string& func_name, aby3::i64Matrix& test, aby3::i64Matrix& res);
template <aby3::Decimal D>
bool check_result(const std::string& func_name, aby3::f64Matrix<D>& test, aby3::f64Matrix<D>& res){
    return check_result(func_name, test.i64Cast(), res.i64Cast());
}