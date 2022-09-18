#include "RTRTest.h"
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>

#include "BuildingBlocks.h"
#include "CipherIndex.h"


using namespace oc;
using namespace aby3;
using namespace std;

// sint and sbit matrix multiplication test.
int test_mul(CLP &cmd){
  IOService ios;
  u64 rows = 10, cols = 1;
  i64Matrix plainInt(rows, cols);
  i64Matrix plainBits(rows, cols);

  for(int i=0; i<rows; i++){
    plainInt(i, 0) = i;
    plainBits(i, 0) = (i%2 == 0) * -1;
  }

  vector<thread> thrds;
  for (int i = 0; i < 3; i++) {
    thrds.emplace_back(thread([i, plainInt, plainBits, &ios]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup((u64)i, ios, enc, eval, runtime); // same setup function in the tutorial.
      si64Matrix ea(plainInt.rows(), plainInt.cols()), eb(plainBits.rows(), plainBits.cols()), res(plainBits.rows(), plainBits.cols());

      i64Matrix zeros(plainInt.rows(), plainInt.cols());
      for(int j=0; j<plainInt.rows(); j++) zeros(j, 0) = 0;
      si64Matrix izeros(zeros.rows(), zeros.cols());

      if (i == 0) {
        enc.localIntMatrix(runtime, plainInt, ea).get();
        enc.localIntMatrix(runtime, plainBits, eb).get();
      } else {
        enc.remoteIntMatrix(runtime, ea).get();
        enc.remoteIntMatrix(runtime, eb).get();
      }
      sbMatrix bitsM;
      Sh3Task task = runtime.noDependencies();
      fetch_msb(i, eb, bitsM, eval, runtime, task);

      // eval.asyncMul(runtime, ea, bitsM, res).get();
      cipher_mul_seq(i, ea, bitsM, res, eval, enc, runtime);
      // cout << "after mul seq" << endl;
      i64Matrix pres;
      enc.revealAll(runtime, res, pres).get();
      if(i == 0){
        cout << "Expected res = [0, 0, 2, 0, 4, 0, 6, ...]" << endl;
        cout << "Real compute = ";
        for(int j=0; j<ea.rows(); j++) cout << pres(j, 0) << " ";
        cout << "\n" << endl;
      }
      // cout << "\n" << endl;
    }));
  }
  for (auto &t : thrds)
    t.join();
  return 0;
}

// // sfixed greater-than test.
int test_gt(CLP &cmd){
  IOService ios;
  u64 rows = 10, cols = 1;
  f64Matrix<D8> plainA(rows, cols);
  f64Matrix<D8> plainB(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      plainA(j, i) = (double) j;
      plainB(j, i) = (double) j + 1;
    }
  }

  vector<thread> thrds;
  for(int i=0; i<3; i++){
    thrds.emplace_back(thread([i, plainA, plainB, &ios]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup((u64)i, ios, enc, eval, runtime); // same setup function in the tutorial.

      sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
      sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());

      if(i == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
      }

      sbMatrix lt_res;
      cipher_gt(i, sharedB, sharedA, lt_res, eval, runtime);
      
      // si64Matrix larger_res;
      sf64Matrix<D8> res;
      res.resize(sharedA.rows(), sharedA.cols());
      eval.asyncMul(runtime, sharedA.i64Cast(), lt_res, res.i64Cast()).get();
      
      f64Matrix<D8> pres;
      enc.revealAll(runtime, res, pres).get();
      
      if(i == 0){
        cout << "Expected res = [0.0, 1.0, 2.0, ..., 9.0]" << endl;
        cout << "Real compute = ";
        for(int j=0; j<sharedA.rows(); j++) cout << pres(j, 0) << ", ";
        cout << "\n" << endl;
      }
    }));
  }

  for (auto &t : thrds)
    t.join();
  return 0;
}

// int test_eq(CLP& cmd){
//   IOService ios;
//   u64 rows = 10, cols = 1;

//   // test for f64 datatype.
//   f64Matrix<D8> plainA(rows, cols);
//   f64Matrix<D8> plainB(rows, cols);
//   f64Matrix<D8> ones(rows, cols);

//   for(u64 j=0; j<rows; j++){
//     for (u64 i = 0; i<cols; i++){
//       plainA(j, i) = (double) j;
//       ones(j, i) = 1;
//       if (j%2 == 0) plainB(j, i) = (double) j + 1;
//       else plainB(j, i) = (double) j;
//     }
//   } 

//   // test for i64 datatype.
//   i64Matrix iplainA(rows, cols);
//   i64Matrix iplainB(rows, cols);
//   i64Matrix iones(rows, cols);

//   for(u64 j=0; j<rows; j++){
//     for (u64 i = 0; i<cols; i++){
//       iplainA(j, i) = j;
//       iones(j, i) = 1;
//       // iplainB(j, i) = j;
//       if (j%2 == 0) iplainB(j, i) = j + 1;
//       else iplainB(j, i) = j;
//     }
//   }

//   vector<thread> thrds;
//   for(u64 i=0; i<3; i++){
//     thrds.emplace_back(thread([i, plainA, plainB, ones, iplainA, iplainB, iones, &ios]() {
//       Sh3Encryptor enc;
//       Sh3Evaluator eval;
//       Sh3Runtime runtime;
//       basic_setup(i, ios, enc, eval, runtime); // same setup function in the tutorial.

//       sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
//       sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());
//       sf64Matrix<D8> sharedOnes(plainA.rows(), plainA.cols());

//       si64Matrix isharedA(plainA.rows(), plainA.cols());
//       si64Matrix isharedB(plainA.rows(), plainA.cols());
//       si64Matrix isharedOnes(plainA.rows(), plainA.cols());

//       if(i == 0){
//         enc.localFixedMatrix(runtime, plainA, sharedA).get();
//         enc.localFixedMatrix(runtime, plainB, sharedB).get();
//         enc.localFixedMatrix(runtime, ones, sharedOnes).get();

//         enc.localIntMatrix(runtime, iplainA, isharedA).get();
//         enc.localIntMatrix(runtime, iplainB, isharedB).get();
//         enc.localIntMatrix(runtime, iplainA, isharedOnes).get();
//       }
//       else{
//         enc.remoteFixedMatrix(runtime, sharedA).get();
//         enc.remoteFixedMatrix(runtime, sharedB).get();
//         enc.remoteFixedMatrix(runtime, sharedOnes).get();

//         enc.remoteIntMatrix(runtime, isharedA).get();
//         enc.remoteIntMatrix(runtime, isharedB).get();
//         enc.remoteIntMatrix(runtime, isharedOnes).get();
//       }

//       sbMatrix eq_res;
//       cipher_eq(i, sharedA, sharedB, eq_res, eval, runtime);
//       sf64Matrix<D8> res;
//       res.resize(sharedA.rows(), sharedA.cols()); 
//       eval.asyncMul(runtime, sharedOnes.i64Cast(), eq_res, res.i64Cast()).get();
//       f64Matrix<D8> pres;
//       enc.revealAll(runtime, res, pres).get();
      
//       if(i == 0){
//         cout << "expected res = [1.0, 0.0, 1.0, ...]" << endl;
//         cout << "equal test = " << pres << endl;
//       }

//       sbMatrix ieq_res;
//       cipher_eq(i, isharedA, isharedB, ieq_res, eval, runtime);

//       si64Matrix ires;
//       ires.resize(isharedA.rows(), isharedA.cols());
//       eval.asyncMul(runtime, isharedOnes, ieq_res, ires).get();
//       i64Matrix ipres;
//       enc.revealAll(runtime, ires, ipres).get();
//       if(i == 0){
//         cout << "test si equal" << endl;
//         cout << "expected res = [1, 0, 1, ...]" << endl;
//         cout << "equal test = " << ipres << endl;
//       }

//     }));
//   }

//   for (auto &t : thrds)
//     t.join();
//   return 0;
// }

int test_argsort(CLP& cmd){
  IOService ios;

  // generate the test data.
  u64 rows = 100, cols = 1;
  f64Matrix<D8> plainTest(rows, cols);
  for(int i=0; i<rows; i++){
    plainTest(i, 0) = i;
  }

  // start the three parties computation.
  vector<thread> thrds;
  for(u64 i=0; i<3; i++){
    thrds.emplace_back(thread([i, &ios, plainTest, rows, cols]() {
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime);

      // generate the cipher test data.
      sf64Matrix<D8> sharedM(rows, cols);
      if(i == 0){
        enc.localFixedMatrix(runtime, plainTest, sharedM).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedM).get();
      }
      
      si64Matrix argres;
      rtr_cipher_argsort(i, sharedM, argres, eval, runtime, enc);
      i64Matrix plain_argres;
      enc.revealAll(runtime, argres, plain_argres).get();
      
      if(i == 0){
        cout << "Expected res = [0 1 2 3 4 ...]" << endl;
        // cout << plain_argres << endl;
        cout << "Real compute = ";
        for(int j=0; j<7; j++){
          cout << plain_argres(j, 0) << " ";
        }
        cout << "\n" << endl;
      }
    }));
  }
  for (auto &t : thrds)
    t.join();
  return 0;
}

int basic_performance(CLP& cmd, int n, int repeats, map<string, vector<double>>& dict){
  IOService ios;
  u64 rows = n, cols = 1;
  f64Matrix<D8> plainA(rows, cols);
  f64Matrix<D8> plainB(rows, cols);

  for(u64 j=0; j<rows; j++){
    for (u64 i = 0; i<cols; i++){
      plainA(j, i) = (double) j;
      plainB(j, i) = (double) j + 1;
    }
  }

  vector<double> time_mul_array(3);
  vector<double> time_gt_array(3);
  vector<double> time_add_array(3);

  vector<thread> thrds;
  for(u64 i=0; i<3; i++){
    thrds.emplace_back(thread([i, plainA, plainB, repeats, &time_mul_array, &time_gt_array, &time_add_array, &ios]() {
      
      // setup the environment.
      clock_t start, end;
      Sh3Encryptor enc;
      Sh3Evaluator eval;
      Sh3Runtime runtime;
      basic_setup(i, ios, enc, eval, runtime);

      // generate the test matrix.
      sf64Matrix<D8> sharedA(plainA.rows(), plainA.cols());
      sf64Matrix<D8> sharedB(plainB.rows(), plainB.cols());

      if(i == 0){
        enc.localFixedMatrix(runtime, plainA, sharedA).get();
        enc.localFixedMatrix(runtime, plainB, sharedB).get();
      }
      else{
        enc.remoteFixedMatrix(runtime, sharedA).get();
        enc.remoteFixedMatrix(runtime, sharedB).get();
      }

      // test the performance of basic ops.
      sbMatrix lt_res;
      start = clock();
      for(int k=0; k<repeats; k++)
        cipher_gt(i, sharedB, sharedA, lt_res, eval, runtime);
      end = clock();

      double time_gt = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
      time_gt_array[i] = time_gt;
      // time_gt_array[i] = 0;
      // // cout << "time_gt_array: " << time_gt_array[i] << endl;

      sf64Matrix<D8> mul_res(plainA.rows(), plainA.cols());
      start = clock();
      for(int k=0; k<repeats; k++)
        // cipher_mul(i, sharedA, sharedB, mul_res, eval, runtime);
        eval.asyncMul(runtime, sharedA, sharedB, mul_res).get();
      end = clock();

      double time_mul = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
      time_mul_array[i] = time_mul;
      // time_mul_array[i] = 0;


      sf64Matrix<D8> add_res(plainA.rows(), plainA.cols());
      start = clock();
      for(int k=0; k<repeats; k++)
        add_res = sharedA + sharedB;
      end = clock();

      double time_add = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
      time_add_array[i] = time_add;

    }));
  }

  for (auto &t : thrds)
    t.join();
  
  // map<string, double(*)[3]> dict;
  dict["mul"] = time_mul_array;
  dict["gt"] = time_gt_array;
  dict["add"] = time_add_array;


  return 0;
}