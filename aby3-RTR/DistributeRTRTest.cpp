#include "DistributeRTRTest.h"

#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>

#include "BuildingBlocks.h"

using namespace oc;
using namespace aby3;
using namespace std;

int dis_test_mul(CLP& cmd){

    // set the role for this process.
    int role = -1;
    if(cmd.isSet("role")){
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if(role == -1){
        throw std::runtime_error(LOCATION);
    }

    // setup the corresponding communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    distribute_setup((u64)role, ios, enc, eval, runtime);

    // run test functions
    // 1. set the test data.
    u64 rows = 10, cols = 1;
    i64Matrix plainInt(rows, cols);
    i64Matrix plainBits(rows, cols);

    for(int i=0; i<rows; i++){
        plainInt(i, 0) = i;
        plainBits(i, 0) = (i%2 == 0) * -1;
    }
    si64Matrix ea(plainInt.rows(), plainInt.cols()), eb(plainBits.rows(), plainBits.cols()), res(plainBits.rows(), plainBits.cols());

    i64Matrix zeros(plainInt.rows(), plainInt.cols());
    for(int j=0; j<plainInt.rows(); j++) zeros(j, 0) = 0;
    si64Matrix izeros(zeros.rows(), zeros.cols());

    if (role == 0) {
        enc.localIntMatrix(runtime, plainInt, ea).get();
        enc.localIntMatrix(runtime, plainBits, eb).get();
    } else {
        enc.remoteIntMatrix(runtime, ea).get();
        enc.remoteIntMatrix(runtime, eb).get();
    }
    sbMatrix bitsM;

    // 2. run function.
    Sh3Task task = runtime.noDependencies();
    fetch_msb(role, eb, bitsM, eval, runtime, task);

    // eval.asyncMul(runtime, ea, bitsM, res).get();
    cipher_mul_seq(role, ea, bitsM, res, eval, enc, runtime);
    // cout << "after mul seq" << endl;
    i64Matrix pres;
    enc.revealAll(runtime, res, pres).get();
    if(role == 0){
        cout << "Expected res = [0, 0, 2, 0, 4, 0, 6, ...]" << endl;
        cout << "Real compute = ";
        for(int j=0; j<ea.rows(); j++) cout << pres(j, 0) << " ";
        cout << "\n" << endl;
    }

    return 0;
}


int dis_basic_performance(CLP& cmd, int n, int repeats, map<std::string, double>& dict){
    // set the role for this process.
    int role = -1;
    if(cmd.isSet("role")){
        auto keys = cmd.getMany<int>("role");
        role = keys[0];
    }
    if(role == -1){
        throw std::runtime_error(LOCATION);
    }

    // setup the corresponding communications.
    IOService ios;
    Sh3Encryptor enc;
    Sh3Evaluator eval;
    Sh3Runtime runtime;
    distribute_setup((u64)role, ios, enc, eval, runtime);

    // set the test data
    u64 rows = n, cols = 1;
    f64Matrix<D8> plainA(rows, cols);
    f64Matrix<D8> plainB(rows, cols);
    i64Matrix iplainA(rows, cols);
    i64Matrix iplainB(rows, cols);

    for(u64 j=0; j<rows; j++){
        for (u64 i = 0; i<cols; i++){
            plainA(j, i) = (double) j;
            plainB(j, i) = (double) j + 1;
            iplainA(j, i) = j;
            iplainB(j, i) = j + 1;
        }
    }

    clock_t start, end;
    double time_mul, time_gt, time_add, time_imul, time_ibmul; // record the total time for each ops in ms.

    sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
    sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());
    si64Matrix isharedA(iplainA.rows(), iplainA.cols());
    si64Matrix isharedB(iplainB.rows(), iplainB.cols());

    if(role == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
        enc.localIntMatrix(runtime, iplainA, isharedA).get();
        enc.localIntMatrix(runtime, iplainB, isharedB).get();
    }
    else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
        enc.remoteIntMatrix(runtime, isharedA).get();
        enc.remoteIntMatrix(runtime, isharedB).get();
    }

    // test the performance of basic ops.
    // 1. sfixed multiplication.
    sf64Matrix<D8> mul_res(plainA.rows(), plainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        eval.asyncMul(runtime, sharedA, sharedB, mul_res).get();
    end = clock();
    time_mul = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);

    // 2. sfixed addition
    sf64Matrix<D8> add_res(plainA.rows(), plainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        add_res = sharedA + sharedB;
    end = clock();
    time_add = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);

    // 3. sfixed greater-then
    sbMatrix lt_res;
    start = clock();
    for(int k=0; k<repeats; k++)
        cipher_gt(role, sharedB, sharedA, lt_res, eval, runtime);
    end = clock();
    time_gt = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);

    // 4. sint multiplication.
    si64Matrix imul_res(iplainA.rows(), iplainB.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        eval.asyncMul(runtime, isharedA, isharedB, imul_res).get();
    end = clock();
    time_imul = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);

    // 5. sb & si multiplication.
    si64Matrix ibmul_res(iplainA.rows(), iplainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        // cipher_mul_seq(role, isharedA, lt_res, ibmul_res, eval, enc, runtime);
        eval.asyncMul(runtime, isharedA, lt_res, ibmul_res).get();
    end = clock();
    time_ibmul = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);

    dict["mul"] = time_mul;
    dict["gt"] = time_gt;
    dict["add"] = time_add;
    dict["imul"] = time_imul;
    dict["ibmul"] = time_ibmul;

    return 0;
}

