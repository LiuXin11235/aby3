#include "PTRProfile.h"
#include <chrono>
#include <thread>
#include "./Pair_then_Reduce/include/datatype.h"
#include "BuildingBlocks.h"
#include "PTRFunction.h"

int profile_cipher_index(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap){

  // set the config info.
  clock_t start, end;
  // set the log file.
  std::string LOG_FOLDER = "/root/aby3/Record/Prof_index/";
  if (cmd.isSet("logFolder")) {
    auto keys = cmd.getMany<std::string>("logFolder");
    LOG_FOLDER = keys[0];
  }
  // cout << "LOG_FOLDER: " << LOG_FOLDER << endl;
  std::string logging_file = LOG_FOLDER + "probe.log";
  std::string profiler_file = LOG_FOLDER + "probe.res";

  // environment setup.
  int role = -1;
  int repeats_ = 100;
  if (cmd.isSet("role")) {
    auto keys = cmd.getMany<int>("role");
    role = keys[0];
  }
  if (role == -1) {
    throw std::runtime_error(LOCATION);
  }
  if (cmd.isSet("repeats")) {
    auto keys = cmd.getMany<int>("repeats");
    repeats_ = keys[0];
  }

  // setup communications.
  IOService ios; Sh3Encryptor enc; Sh3Evaluator eval; Sh3Runtime runtime;
  distribute_setup((u64)role, ios, enc, eval, runtime);


  // task setup.
  auto mpiPtrTask =
      new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(
          1, vector_size_start, role, enc, runtime, eval);
  mpiPtrTask->circuit_construct({(size_t)vector_size_start}, {(size_t)n});

  // data construct.
  m = vector_size_start;
  aby3::si64 dval; dval.mData[0] = 0; dval.mData[1] = 0;

  size_t m_start = mpiPtrTask->m_start; size_t m_end = mpiPtrTask->m_end;
  size_t partial_len = m_end - m_start + 1;

  vector<si64> res(m); vector_generation(role, enc, runtime, res);
  vector<si64> vecM(partial_len); vector_generation(role, enc, runtime, vecM);
  vector<si64> vecIndex(m); vector_generation(role, enc, runtime, vecIndex);
  vector<int> range_index(partial_len); vector_generation(role, enc, runtime, range_index);
  
  // construct the fake data
  mpiPtrTask->set_selective_value(vecM.data(), 0);
  FakeArray<aby3::si64> dataX =
      fake_repeat(vecIndex.data(), mpiPtrTask->shapeX, mpiPtrTask->m, 0);
  FakeArray<int> dataY = fake_repeat(
      range_index.data(), mpiPtrTask->shapeY, mpiPtrTask->n, 1, mpiPtrTask->m_start,
      mpiPtrTask->m_end - mpiPtrTask->m_start + 1);

  // begin the profiler: 0) set the initial
  start = clock();
  for (int i = 0; i < repeats_; i++) mpiPtrTask->subTask->circuit_profile(dataX, dataY, mpiPtrTask->selectV);
  end = clock();

  // time synchronized.
  double start_time = double((end - start) * 1000) / (repeats_ * CLOCKS_PER_SEC);
  start_time = synchronized_time(role, start_time, runtime);

  double double_time = -1;
  size_t vector_size = vector_size_start;
  double last_ratio = start_time / vector_size; double ratio = last_ratio / 2;

  // first: 1) exp
  while(ratio < last_ratio){
    vector_size *= 2;
    auto testMpiTask =
        new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(
            1, vector_size, role, enc, runtime, eval);
    testMpiTask->circuit_construct({vector_size}, {n});
    m = vector_size;

    // data construct
    size_t m_start_ = testMpiTask->m_start;
    size_t m_end_ = testMpiTask->m_end;
    size_t partial_len_ = m_end_ - m_start_ + 1;

    vector<si64> res_(m); vector_generation(role, enc, runtime, res_);
    vector<si64> vecM_(partial_len_); vector_generation(role, enc, runtime, vecM_);
    vector<si64> vecIndex_(m); vector_generation(role, enc, runtime, vecIndex_);
    vector<int> range_index_(partial_len_); vector_generation(role, enc, runtime, range_index_);

    FakeArray<aby3::si64> dataX_ =
        fake_repeat(vecIndex_.data(), testMpiTask->shapeX, testMpiTask->m, 0);
    FakeArray<int> dataY_ = fake_repeat(
        range_index_.data(), testMpiTask->shapeY, testMpiTask->n, 1,
        testMpiTask->m_start, testMpiTask->m_end - testMpiTask->m_start + 1);

    testMpiTask->set_default_value(dval);
    testMpiTask->set_selective_value(vecM_.data(), 0);

#ifdef LOGING
    write_log(logging_file,
              "!!!! staring with vector: " + to_string(vector_size));
    this_thread::sleep_for(chrono::seconds(3));
#endif

    // adjusting the repeat time.
    if (vector_size > 5000 && vector_size < 50000) repeats_ = 1000;
    if (vector_size > 50000) repeats_ = 50;

    start = clock();
    for (int i = 0; i < repeats_; i++) {
      testMpiTask->subTask->circuit_profile(dataX_, dataY_,
                                            testMpiTask->selectV);
    }
    end = clock();
    double_time = (end-start)*1000 / (repeats_ * CLOCKS_PER_SEC);
    double_time = synchronized_time(role, double_time, runtime);

    last_ratio = ratio;
    ratio = double_time / vector_size;

#ifdef LOGING
    write_log(logging_file, "time = " + to_string(double_time) +
                                " | vector_size = " + to_string(vector_size) +
                                "| ratio = " +
                                to_string(double_time / vector_size));
    write_log(logging_file,
              "time diff = " + to_string(double_time - start_time));
#endif
  }

  // second: 2) binary search
  size_t vector_start = vector_size / 2;
  size_t vector_end = vector_size;

  while ((vector_end - vector_start) > gap) {
    vector_size = (size_t)(vector_start + vector_end) / 2;

    // measure time
    auto testMpiTask =
        new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(
            1, vector_size, role, enc, runtime, eval);
    testMpiTask->circuit_construct({vector_size}, {n});
    m = vector_size;

    // data construct
    size_t m_start_ = testMpiTask->m_start;
    size_t m_end_ = testMpiTask->m_end;
    size_t partial_len_ = m_end_ - m_start_ + 1;

    vector<si64> res_(m); vector_generation(role, enc, runtime, res_);
    vector<si64> vecM_(partial_len_); vector_generation(role, enc, runtime, vecM_);
    vector<si64> vecIndex_(m); vector_generation(role, enc, runtime, vecIndex_);
    vector<int> range_index_(partial_len_); vector_generation(role, enc, runtime, range_index_);

    FakeArray<aby3::si64> dataX_ =
        fake_repeat(vecIndex_.data(), testMpiTask->shapeX, testMpiTask->m, 0);
    FakeArray<int> dataY_ = fake_repeat(
        range_index_.data(), testMpiTask->shapeY, testMpiTask->n, 1,
        testMpiTask->m_start, testMpiTask->m_end - testMpiTask->m_start + 1);

    testMpiTask->set_default_value(dval);
    testMpiTask->set_selective_value(vecM_.data(), 0);

    start = clock();
    for (int i = 0; i < repeats_; i++) {
      // cout << "in test: " << i << endl;
      testMpiTask->subTask->circuit_profile(dataX_, dataY_,
                                            testMpiTask->selectV);
    }
    end = clock();
    double_time = double((end - start)*1000)/(repeats_ * CLOCKS_PER_SEC);
    synchronized_time(role, double_time, runtime);

    last_ratio = ratio;
    ratio = double_time / vector_size;

#ifdef LOGING
    write_log(logging_file, "time = " + to_string(double_time) +
                                " | vector_size = " + to_string(vector_size) +
                                "| ratio = " +
                                to_string(double_time / vector_size));
#endif

    if (ratio > last_ratio) vector_end = vector_size;
    else vector_start = vector_size;
  }

  size_t optimal_b = (size_t)((vector_start + vector_end) / 2);

  // save the optimalB to config file.
  std::ofstream resfs(profiler_file, std::ios_base::app);
  resfs << "optimal_B: " << optimal_b << std::endl;
  resfs.close();
  
  return 0;
}