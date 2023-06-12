#include <chrono>  
#include <thread> 
#include "PTRTest.h"
#include "PTRFunction.h"
#include "BuildingBlocks.h"
#include "./Pair_then_Reduce/include/datatype.h"

using namespace oc;
using namespace aby3;
using namespace std;

// #define DEBUG
#define LOGING

int test_cipher_index_ptr(CLP& cmd, size_t n, size_t m){

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  distribute_setup((u64)role, ios, enc, eval, runtime);

  // generate the test data.
  i64Matrix plainTest(n, 1);
  i64Matrix plainIndex(m, 1);
  for (int i = 0; i < n; i++) {
    plainTest(i, 0) = i;
  }
  // inverse sequence.
  for (int i = 0; i < m; i++) {
    plainIndex(i, 0) = n - 1 - i;
  }

  // generate the cipher test data.
  si64Matrix sharedM(n, 1);
  si64Matrix sharedIndex(m, 1);
  if(role == 0){
      enc.localIntMatrix(runtime, plainTest, sharedM).get();
      enc.localIntMatrix(runtime, plainIndex, sharedIndex).get();
  }
  else{
      enc.remoteIntMatrix(runtime, sharedM).get();
      enc.remoteIntMatrix(runtime, sharedIndex).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecM(n);
  vector<si64> vecIndex(m);
  for(int i=0; i<m; i++) vecIndex[i] = sharedIndex(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<n; i++) vecM[i] = sharedM(i, 0);

  ptr_secret_index(role, vecM, vecIndex, res, eval, runtime, enc);

  // test for correctness:
  si64Matrix aby3res(m, 1);
  for(int i=0; i<m; i++){
    aby3res(i, 0, res[i]);
  }

  i64Matrix plain_res;
  enc.revealAll(runtime, aby3res, plain_res).get();
  if (role == 0) {
    for(int i=0; i<m; i++){
      if(plain_res(i, 0) != plainTest(n - 1 - i)){
        cout << "res != test in index " << i << "\nres[i] = " <<  plain_res(i, 0) << " plainTest = " << plainTest(n - 1 - i) << endl;
      }
    }
    cout << "FINISH!" << endl;
  }
}


int test_cipher_index_ptr_mpi(CLP& cmd, size_t n, size_t m, int task_num, int opt_B){

  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  // cout << "in this function" << endl;

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_index/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(m) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);


  start = clock();
  // construct task
  auto mpiPtrTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(task_num, opt_B, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({m}, {n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;

  // directly generate the partial data
  size_t partial_len = m_end - m_start + 1;
  i64Matrix plainTest(partial_len, 1);
  for (int i = 0; i < partial_len; i++) {
    plainTest(i, 0) = i + m_start;
  }
  int* range_index = new int[partial_len];
  for(int i=m_start; i<m_end + 1; i++) range_index[i-m_start] = i;


  i64Matrix plainIndex(m, 1);
  // inverse sequence.
  for (int i = 0; i < m; i++) {
    plainIndex(i, 0) = n - 1 - i;
  }

  // generate the cipher test data.
  si64Matrix sharedM(partial_len, 1);
  si64Matrix sharedIndex(m, 1);
  if(role == 0){
      enc.localIntMatrix(runtime, plainTest, sharedM).get();
      enc.localIntMatrix(runtime, plainIndex, sharedIndex).get();
  }
  else{
      enc.remoteIntMatrix(runtime, sharedM).get();
      enc.remoteIntMatrix(runtime, sharedIndex).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecM(partial_len);
  vector<si64> vecIndex(m);
  for(int i=0; i<m; i++) vecIndex[i] = sharedIndex(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecM[i] = sharedM(i, 0);

  // set the data related part.
  mpiPtrTask->set_selective_value(vecM.data(), 0);
  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  if(rank == 0){
    cout << "binning circuit evaluate..." << endl;
  }
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(vecIndex.data(), range_index, vecM.data(), res.data());

  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrTask->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrTask->time_combine << "\n" << std::endl;
    ofs.close();
  }
}


int test_cipher_select_ptr_mpi(CLP& cmd, size_t n, size_t m, int task_num, int opt_B){

  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  // cout << rank << endl;

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_select/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(m) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);


  start = clock();
  // construct task
  auto mpiPtrTask = new MPISecretIndex<int, aby3::si64, aby3::si64, aby3::si64, SubIndex>(task_num, opt_B, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({m}, {n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;

  // directly generate the partial data
  size_t partial_len = m_end - m_start + 1;
  i64Matrix plainTest(partial_len, 1);
  for (int i = 0; i < partial_len; i++) {
    plainTest(i, 0) = i + m_start;
  }
  int* range_index = new int[m];
  for(int i=0; i<m; i++) range_index[i] = i;


  i64Matrix plainIndex(partial_len, 1);
  // inverse sequence.
  for (int i = 0; i < partial_len; i++) {
    plainIndex(i, 0) = n - 1 - (i + m_start);
  }

  // generate the cipher test data.
  si64Matrix sharedM(partial_len, 1);
  si64Matrix sharedIndex(partial_len, 1);
  if(role == 0){
      enc.localIntMatrix(runtime, plainTest, sharedM).get();
      enc.localIntMatrix(runtime, plainIndex, sharedIndex).get();
  }
  else{
      enc.remoteIntMatrix(runtime, sharedM).get();
      enc.remoteIntMatrix(runtime, sharedIndex).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecM(partial_len);
  vector<si64> vecIndex(partial_len);
  for(int i=0; i<partial_len; i++) vecIndex[i] = sharedIndex(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecM[i] = sharedM(i, 0);

  // set the data related part.
  mpiPtrTask->set_selective_value(vecM.data(), 0);
  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(range_index, vecIndex.data(), vecM.data(), res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrTask->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrTask->time_combine << "\n" << std::endl;
    ofs.close();
  }
}


int test_cipher_rank_ptr_mpi(oc::CLP& cmd, size_t n, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  // cout << rank << endl;

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_rank/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(n) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);


  start = clock();
  // construct task
  auto mpiPtrTask = new MPIRank<aby3::si64, aby3::si64, aby3::si64, aby3::si64, SubRank>(task_num, opt_B, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({(size_t) n}, {(size_t) n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;

  // generate raw test data
  i64Matrix testData(n, 1);
  for(int i=0; i<n; i++) testData(i, 0) = i;

  size_t partial_len = m_end - m_start + 1;
  i64Matrix partialData(partial_len, 1);
  for(int i=0; i<partial_len; i++) partialData(i, 0) = (i + m_start);

  // generate cipher test data
  si64Matrix sharedM(n, 1);
  si64Matrix partialM(partial_len, 1);
  if(role == 0){
    enc.localIntMatrix(runtime, testData, sharedM).get();
    enc.localIntMatrix(runtime, partialData, partialM).get();
  }
  else{
    enc.remoteIntMatrix(runtime, sharedM).get();
    enc.remoteIntMatrix(runtime, partialM).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, n);

  vector<si64> res(n);
  vector<si64> vecM(n);
  vector<si64> vecPartialM(partial_len);
  for(int i=0; i<n; i++) vecM[i] = sharedM(i, 0);
  for(int i=0; i<n; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecPartialM[i] = partialM(i, 0);

  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(vecM.data(), vecPartialM.data(), nullptr, res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);


  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrTask->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrTask->time_combine << "\n" << std::endl;
    ofs.close();
  }
}


int test_cipher_sort_ptr_mpi(oc::CLP& cmd, size_t n, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  // cout << rank << endl;

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_sort/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(n) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  // construct task
  auto mpiPtrRank = new MPIRank<aby3::si64, aby3::si64, aby3::si64, aby3::si64, SubRank>(task_num, opt_B, role, enc, runtime, eval);
  auto mpiPtrIndex = new MPISecretIndex<int, aby3::si64, aby3::si64, aby3::si64, SubIndex>(task_num, opt_B, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrRank->set_default_value(dval);
  mpiPtrRank->circuit_construct({(size_t) n}, {(size_t) n});
  mpiPtrIndex->set_default_value(dval);
  mpiPtrIndex->circuit_construct({(size_t) n}, {(size_t) n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrRank->m_start;
  size_t m_end = mpiPtrRank->m_end;

  // generate raw test data
  i64Matrix testData(n, 1);
  for(int i=0; i<n; i++) testData(i, 0) = i;

  size_t partial_len = m_end - m_start + 1;
  i64Matrix partialData(partial_len, 1);
  for(int i=0; i<partial_len; i++) partialData(i, 0) = (i + m_start);

  // generate cipher test data
  si64Matrix sharedM(n, 1);
  si64Matrix partialM(partial_len, 1);
  if(role == 0){
    enc.localIntMatrix(runtime, testData, sharedM).get();
    enc.localIntMatrix(runtime, partialData, partialM).get();
  }
  else{
    enc.remoteIntMatrix(runtime, sharedM).get();
    enc.remoteIntMatrix(runtime, partialM).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, n);

  vector<si64> res(n);
  vector<si64> vecM(n);
  vector<si64> vecPartialM(partial_len);
  for(int i=0; i<n; i++) vecM[i] = sharedM(i, 0);
  for(int i=0; i<n; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecPartialM[i] = partialM(i, 0);

  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  // evaluate the task.
  mpiPtrRank->circuit_evaluate(vecM.data(), vecPartialM.data(), nullptr, res.data());

  // TODO:
  // mpiPtrIndex->data();


  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrRank->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrRank->time_combine << "\n" << std::endl;
    ofs.close();
  }
}


int test_cipher_search_ptr_mpi(oc::CLP& cmd, size_t n, size_t m, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  // cout << rank << endl;

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_search/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(m) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);


  start = clock();
  // construct task
  auto mpiPtrTask = new MPISearch<aby3::si64, aby3::si64, aby3::si64, aby3::si64, SubSearch>(task_num, opt_B, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({m}, {n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;

  // directly generate the partial data
  size_t partial_len = m_end - m_start + 1;
  i64Matrix plainTest(partial_len, 1);
  for (int i = 0; i < partial_len; i++) {
    plainTest(i, 0) = i + m_start;
  }

  // non-sense diffValue, just for test.
  i64Matrix diffValue(partial_len, 1);
  for(int i=0; i<partial_len; i++){
    diffValue(i, 0) = i + m_start;
  }


  i64Matrix searchKey(m, 1);
  // inverse sequence.
  for (int i = 0; i < m; i++) {
    searchKey(i, 0) = n - 1 - i;
  }

  // generate the cipher test data.
  si64Matrix secretDiff(partial_len, 1);
  si64Matrix secretSpace(partial_len, 1);
  si64Matrix secretKey(m, 1);
  if(role == 0){
      enc.localIntMatrix(runtime, plainTest, secretSpace).get();
      enc.localIntMatrix(runtime, diffValue, secretDiff).get();
      enc.localIntMatrix(runtime, searchKey, secretKey).get();
  }
  else{
      enc.remoteIntMatrix(runtime, secretSpace).get();
      enc.remoteIntMatrix(runtime, secretDiff).get();
      enc.remoteIntMatrix(runtime, secretKey).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecSpace(partial_len);
  vector<si64> vecDiff(partial_len);
  vector<si64> vecKey(m);
  for(int i=0; i<m; i++) vecKey[i] = secretKey(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++){
    vecDiff[i] = secretDiff(i, 0);
    vecSpace[i] = secretSpace(i, 0);
  }

  // set the data related part.
  mpiPtrTask->set_selective_value(vecDiff.data(), 0);
  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(vecKey.data(), vecSpace.data(), vecDiff.data(), res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrTask->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrTask->time_combine << "\n" << std::endl;
    ofs.close();
  }

  return 0;
}


int test_cipher_average_ptr_mpi(oc::CLP& cmd, size_t n, size_t m, int task_num, int opt_B){
  
  //1.  task setup.
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  clock_t start, end;

  // 1) set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_average/";
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(m) + "-M=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-OPT_B=" + std::to_string(opt_B) + "-" + std::to_string(rank);

  int role = -1;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }

  start = clock();
  // 2) setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);

  // 3) initial task.
  start = clock();
  auto mpiPtrTask = new MPIAverage<int, aby3::si64, aby3::si64, aby3::si64, SubAvg>(task_num, opt_B, role, enc, runtime, eval);
  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({m}, {n});
  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  // 2. data generation.
  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;
  size_t partial_len = m_end - m_start + 1;
  si64Matrix data(partial_len, 1);
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  init_ones(role, enc, runtime, data, partial_len);  

  vector<si64> vecM(partial_len);
  vector<int> x(1);
  vector<si64> res(m);
  x[0] = 0;
  for(int i=0; i<partial_len; i++) vecM[i] = data(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  // 3. task evaluate.
  start = clock();
  mpiPtrTask->circuit_evaluate(x.data(), vecM.data(), nullptr, res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << "subTask: " << std::setprecision(5) << mpiPtrTask->time_subTask << "\ntime_combine: " << std::setprecision(5) << mpiPtrTask->time_combine << "\n" << std::endl;
    ofs.close();
  }
}


int profile_index(oc::CLP& cmd, size_t n, size_t m, int vector_size, int task_num){
    // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  clock_t start, end;

  // set the log file.
  std::string LOG_FOLDER="/root/aby3/Record/Record_profile/";
  if(cmd.isSet("logFolder")){
      auto keys = cmd.getMany<std::string>("logFolder");
      LOG_FOLDER = keys[0];
  }
  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(m) + "-TASKS=" + std::to_string(task_num) + "-Vec=" + std::to_string(vector_size) + "-" + std::to_string(rank);

  int role = -1;
  int repeats_ = 1000;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }
  if(cmd.isSet("repeats")){
    auto keys = cmd.getMany<int>("repeats");
    repeats_ = keys[0];
  }

  start = clock();
  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);


  start = clock();
  // construct task
  auto mpiPtrTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(task_num, vector_size, role, enc, runtime, eval);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  mpiPtrTask->set_default_value(dval);
  mpiPtrTask->circuit_construct({m}, {n});

  end = clock(); // time for task init.
  double time_task_init = double((end - start)*1000)/(CLOCKS_PER_SEC);

  start = clock();
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;

  // directly generate the partial data
  size_t partial_len = m_end - m_start + 1;
  i64Matrix plainTest(partial_len, 1);
  for (int i = 0; i < partial_len; i++) {
    plainTest(i, 0) = i + m_start;
  }
  int* range_index = new int[partial_len];
  for(int i=m_start; i<m_end + 1; i++) range_index[i-m_start] = i;


  i64Matrix plainIndex(m, 1);
  // inverse sequence.
  for (int i = 0; i < m; i++) {
    plainIndex(i, 0) = n - 1 - i;
  }

  // generate the cipher test data.
  si64Matrix sharedM(partial_len, 1);
  si64Matrix sharedIndex(m, 1);
  if(role == 0){
      enc.localIntMatrix(runtime, plainTest, sharedM).get();
      enc.localIntMatrix(runtime, plainIndex, sharedIndex).get();
  }
  else{
      enc.remoteIntMatrix(runtime, sharedM).get();
      enc.remoteIntMatrix(runtime, sharedIndex).get();
  }
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecM(partial_len);
  vector<si64> vecIndex(m);
  for(int i=0; i<m; i++) vecIndex[i] = sharedIndex(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecM[i] = sharedM(i, 0);

  // set the data related part.
  mpiPtrTask->set_selective_value(vecM.data(), 0);
  end = clock();
  double time_task_prep = double((end - start)*1000)/(CLOCKS_PER_SEC);

  // construct the fake data
  FakeArray<aby3::si64> dataX = fake_repeat(vecIndex.data(), mpiPtrTask->shapeX, mpiPtrTask->m, 0);
  FakeArray<int> dataY = fake_repeat(range_index, mpiPtrTask->shapeY, mpiPtrTask->n, 1, mpiPtrTask->m_start, mpiPtrTask->m_end - mpiPtrTask->m_start + 1);

  start = clock();
  // evaluate the task.
  for(int i=0; i<repeats_; i++){
    mpiPtrTask->subTask->circuit_profile(dataX, dataY, mpiPtrTask->selectV);
  }
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC*repeats_);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_profile: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }
}


int probe_profile_index(oc::CLP& cmd, size_t n, size_t m, int vector_size_start, double epsilon, size_t gap){

  // used to profile single task with limited bandwidth and see its performence.
  clock_t start, end;
  // cout << "in prob profile index " << endl;
  // set the log file.
  std::string LOG_FOLDER="/root/aby3/Record/Record_probe/";
  if(cmd.isSet("logFolder")){
      auto keys = cmd.getMany<std::string>("logFolder");
      LOG_FOLDER = keys[0];
  }
  // cout << "LOG_FOLDER: " << LOG_FOLDER << endl;
  std::string logging_file = LOG_FOLDER + "probe.log";
  std::string profiler_file = LOG_FOLDER + "probe.res";

  int role = -1;
  int repeats_ = 100;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }
  if(cmd.isSet("repeats")){
    auto keys = cmd.getMany<int>("repeats");
    repeats_ = keys[0];
  }

  // setup communications.
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  distribute_setup((u64)role, ios, enc, eval, runtime);

  aby3::si64 dval;
  dval.mData[0] = 0, dval.mData[1] = 0;
  auto mpiPtrTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(1, vector_size_start, role, enc, runtime, eval);
  mpiPtrTask->circuit_construct({(size_t)vector_size_start}, {(size_t)n});
  m = vector_size_start;
  // data construct
  size_t m_start = mpiPtrTask->m_start;
  size_t m_end = mpiPtrTask->m_end;
  size_t partial_len = m_end - m_start + 1;
  si64Matrix sharedM(partial_len, 1);
  si64Matrix sharedIndex(m, 1);
  init_ones(role, enc, runtime, sharedM, partial_len);  
  init_ones(role, enc, runtime, sharedIndex, m);  
  int* range_index = new int[partial_len];
  for(int i=m_start; i<m_end + 1; i++) range_index[i-m_start] = i;
  si64Matrix init_res;
  init_zeros(role, enc, runtime, init_res, m);
  
  vector<si64> res(m);
  vector<si64> vecM(partial_len);
  vector<si64> vecIndex(m);
  for(int i=0; i<m; i++) vecIndex[i] = sharedIndex(i, 0);
  for(int i=0; i<m; i++) res[i] = init_res(i, 0);
  for(int i=0; i<partial_len; i++) vecM[i] = sharedM(i, 0);

  // construct the fake data
  mpiPtrTask->set_selective_value(vecM.data(), 0);
  FakeArray<aby3::si64> dataX = fake_repeat(vecIndex.data(), mpiPtrTask->shapeX, mpiPtrTask->m, 0);
  FakeArray<int> dataY = fake_repeat(range_index, mpiPtrTask->shapeY, mpiPtrTask->n, 1, mpiPtrTask->m_start, mpiPtrTask->m_end - mpiPtrTask->m_start + 1);

  start = clock();
  for(int i=0; i<repeats_; i++){
    mpiPtrTask->subTask->circuit_profile(dataX, dataY, mpiPtrTask->selectV);
  }
  end = clock();
  double start_time;

  if(role == 0){
    start_time = double((end - start)*1000)/(repeats_ * CLOCKS_PER_SEC);
    runtime.mComm.mNext.asyncSend<double>(start_time);
    runtime.mComm.mPrev.asyncSend<double>(start_time);
  }
  if(role == 1){
    auto tmp = runtime.mComm.mPrev.asyncRecv<double>(&start_time, 1);
    tmp.get();
  }
  if(role == 2){
    auto tmp = runtime.mComm.mNext.asyncRecv<double>(&start_time, 1);
    tmp.get();
  }

  double double_time = -1;
  size_t vector_size = vector_size_start;
  double last_ratio = start_time / vector_size;
  double ratio = last_ratio / 2;

  #ifdef LOGING
  write_log(logging_file, "time = " + to_string(start_time) + " | vector_size = " + to_string(vector_size) + "| ratio = " + to_string(start_time / vector_size));
  #endif

  // exponential probe
  // while(ratio < last_ratio){
  while(true){
    vector_size *= 2;
    auto testMpiTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(1, vector_size, role, enc, runtime, eval);
    testMpiTask->circuit_construct({vector_size}, {n});
    m = vector_size;
    // data construct
    size_t m_start_ = testMpiTask->m_start;
    size_t m_end_ = testMpiTask->m_end;
    size_t partial_len_ = m_end_ - m_start_ + 1;
    si64Matrix sharedM_(partial_len_, 1);
    si64Matrix sharedIndex_(vector_size, 1);
    init_ones(role, enc, runtime, sharedM_, partial_len_);  
    init_ones(role, enc, runtime, sharedIndex_, vector_size);  
    int* range_index_ = new int[partial_len_];
    for(int i=m_start_; i<m_end_ + 1; i++) range_index_[i-m_start_] = i;
    si64Matrix init_res_;
    init_zeros(role, enc, runtime, init_res_, vector_size);
    
    vector<si64> res_(m);
    vector<si64> vecM_(partial_len_);
    vector<si64> vecIndex_(m);
    for(int i=0; i<m; i++) vecIndex_[i] = sharedIndex_(i, 0);
    for(int i=0; i<m; i++) res_[i] = init_res_(i, 0);
    for(int i=0; i<partial_len_; i++) vecM_[i] = sharedM_(i, 0);

    FakeArray<aby3::si64> dataX_ = fake_repeat(vecIndex_.data(), testMpiTask->shapeX, testMpiTask->m, 0);
    FakeArray<int> dataY_ = fake_repeat(range_index_, testMpiTask->shapeY, testMpiTask->n, 1, testMpiTask->m_start, testMpiTask->m_end - testMpiTask->m_start + 1);

    testMpiTask->set_default_value(dval);
    testMpiTask->set_selective_value(vecM_.data(), 0);

    #ifdef LOGING
    write_log(logging_file, "!!!! staring with vector: " + to_string(vector_size));
    this_thread::sleep_for(chrono::seconds(3));
    #endif

    if(vector_size > 5000 && vector_size < 50000){
      repeats_ = 1000;
    }

    if(vector_size > 50000){
      repeats_ = 50;
    }
    
    start = clock();
    for(int i=0; i<repeats_; i++){
      testMpiTask->subTask->circuit_profile(dataX_, dataY_, testMpiTask->selectV);
    }
    end = clock();

    // #ifdef LOGING
    // write_log(logging_file, "time = " + to_string(double_time) + " | vector_size = " + to_string(vector_size));
    // write_log(logging_file, "time diff = " + to_string(double_time - start_time));
    // #endif

    // synchronize with double time...
    if(role == 0){
      double_time = double((end - start)*1000)/(repeats_ * CLOCKS_PER_SEC);
      runtime.mComm.mNext.asyncSend<double>(double_time);
      runtime.mComm.mPrev.asyncSend<double>(double_time);
    }
    if(role == 1){
      auto tmp = runtime.mComm.mPrev.asyncRecv<double>(&double_time, 1);
      tmp.get();
    }
    if(role == 2){
      auto tmp = runtime.mComm.mNext.asyncRecv<double>(&double_time, 1);
      tmp.get();
    }

    last_ratio = ratio;
    ratio = double_time / vector_size;

    #ifdef LOGING
    write_log(logging_file, "time = " + to_string(double_time) + " | vector_size = " + to_string(vector_size) + "| ratio = " + to_string(double_time / vector_size));
    write_log(logging_file, "time diff = " + to_string(double_time - start_time));
    #endif

  }

  // then find B by binary split.
  size_t vector_start = vector_size / 2;
  size_t vector_end = vector_size;

  #ifdef LOGING
  write_log(logging_file, "after exponential probe, interval is at : " + to_string(vector_start) + " to " + to_string(vector_end));
  #endif


  // size_t mid_point;
  while((vector_end - vector_start) > gap){
    vector_size = (size_t) (vector_start + vector_end) / 2;

    // measure time
    auto testMpiTask = new MPISecretIndex<aby3::si64, int, aby3::si64, aby3::si64, SubIndex>(1, vector_size, role, enc, runtime, eval);
    testMpiTask->circuit_construct({vector_size}, {n});
    m = vector_size;

    // data construct
    size_t m_start_ = testMpiTask->m_start;
    size_t m_end_ = testMpiTask->m_end;
    size_t partial_len_ = m_end_ - m_start_ + 1;
    si64Matrix sharedM_(partial_len_, 1);
    si64Matrix sharedIndex_(vector_size, 1);
    init_ones(role, enc, runtime, sharedM_, partial_len_);  
    init_ones(role, enc, runtime, sharedIndex_, vector_size);  
    int* range_index_ = new int[partial_len_];
    for(int i=m_start_; i<m_end_ + 1; i++) range_index_[i-m_start_] = i;
    si64Matrix init_res_(vector_size, 1);
    init_ones(role, enc, runtime, init_res_, vector_size);

    vector<si64> res_(m);
    vector<si64> vecM_(partial_len_);
    vector<si64> vecIndex_(m);
    for(int i=0; i<m; i++) vecIndex_[i] = sharedIndex_(i, 0);
    for(int i=0; i<m; i++) res_[i] = init_res_(i, 0);
    for(int i=0; i<partial_len_; i++) vecM_[i] = sharedM_(i, 0);

    FakeArray<aby3::si64> dataX_ = fake_repeat(vecIndex_.data(), testMpiTask->shapeX, testMpiTask->m, 0);
    FakeArray<int> dataY_ = fake_repeat(range_index_, testMpiTask->shapeY, testMpiTask->n, 1, testMpiTask->m_start, testMpiTask->m_end - testMpiTask->m_start + 1);

    testMpiTask->set_default_value(dval);
    testMpiTask->set_selective_value(vecM_.data(), 0);

    start = clock();
    for(int i=0; i<repeats_; i++){
      // cout << "in test: " << i << endl;
      testMpiTask->subTask->circuit_profile(dataX_, dataY_, testMpiTask->selectV);
    }
    end = clock();
    // double_time = double((end - start)*1000)/(repeats_ * CLOCKS_PER_SEC);

    // synchronize with double time...
    if(role == 0){
      double_time = double((end - start)*1000)/(repeats_ * CLOCKS_PER_SEC);
      runtime.mComm.mNext.asyncSend<double>(double_time);
      runtime.mComm.mPrev.asyncSend<double>(double_time);
    }
    if(role == 1){
      auto tmp = runtime.mComm.mPrev.asyncRecv<double>(&double_time, 1);
      tmp.get();
    }
    if(role == 2){
      auto tmp = runtime.mComm.mNext.asyncRecv<double>(&double_time, 1);
      tmp.get();
    }

    last_ratio = ratio;
    ratio = double_time / vector_size;

    #ifdef LOGING
    write_log(logging_file, "time = " + to_string(double_time) + " | vector_size = " + to_string(vector_size) + "| ratio = " + to_string(double_time / vector_size));
    #endif

    if(ratio > last_ratio) vector_end = vector_size;
    else vector_start = vector_size;
  }

  size_t optimal_b = (size_t) ((vector_start + vector_end) / 2);
  
  // save the optimalB to config file.
  std::ofstream resfs(profiler_file, std::ios_base::app);
  resfs << "optimal_B: " << optimal_b << std::endl;  
  resfs.close();

  return 0;
}


int test_vectorization(oc::CLP& cmd, size_t n_, int task_num){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  clock_t start, end;

  // set the log file.
  static std::string LOG_FOLDER="/root/aby3/Record/Record_vector/";

  int role = -1;
  int repeats_ = 1000;
  if(cmd.isSet("role")){
      auto keys = cmd.getMany<int>("role");
      role = keys[0];
  }
  if(role == -1){
      throw std::runtime_error(LOCATION);
  }
  if(cmd.isSet("repeats")){
    auto keys = cmd.getMany<int>("repeats");
    repeats_ = keys[0];
  }

  // multi-processor deployment setup.
  start = clock();
  IOService ios;
  Sh3Encryptor enc;
  Sh3Evaluator eval;
  Sh3Runtime runtime;
  multi_processor_setup((u64)role, rank, ios, enc, eval, runtime);
  end = clock();
  double time_task_setup = double((end - start)*1000)/(CLOCKS_PER_SEC);

  size_t n = n_;
  int repeats = repeats_;

  std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-" + std::to_string(rank);
  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "begin" << std::endl;  
    ofs.close();
  }

  map<string, double> dict;

  // data preparation.
  start = clock();
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
  end = clock();
  double time_data_prepare = double((end - start)*1000) / CLOCKS_PER_SEC;
  for(int r=0; r<task_num; r++){
    if(rank == r){
        std::ofstream ofs(logging_file, std::ios_base::app);
        ofs << "data_prepare: " << time_data_prepare << std::endl;  
        ofs << "time_task_setup: " << time_task_setup << std::endl;  
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // begin test functions.
  // 1. sfixed multiplication.
  sf64Matrix<D8> mul_res(plainA.rows(), plainA.cols());
  start = clock();
  for(int k=0; k<repeats; k++)
      eval.asyncMul(runtime, sharedA, sharedB, mul_res).get();
  end = clock();
  dict["fxp-mul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "fxp-mul: " << dict["fxp-mul"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "fxp-mul: " << dict["fxp-mul"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // 2. sfixed addition
  sf64Matrix<D8> add_res(plainA.rows(), plainA.cols());
  start = clock();
  for(int k=0; k<repeats; k++)
      add_res = sharedA + sharedB;
  end = clock();
  dict["fxp-add"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "fxp-add: " << dict["fxp-add"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "fxp-add: " << dict["fxp-add"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // 3. sfixed greater-then
  sbMatrix lt_res;
  start = clock();
  for(int k=0; k<repeats; k++)
      cipher_gt(role, sharedB, sharedA, lt_res, eval, runtime);
  end = clock();
  dict["fxp-gt"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "fxp-gt: " << dict["fxp-gt"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "fxp-gt: " << dict["fxp-gt"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  sbMatrix eq_res;
  start = clock();
  for(int k=0; k<repeats; k++)
      cipher_eq(role, sharedB, sharedA, eq_res, eval, runtime);
  end = clock();
  dict["fxp-eq"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "fxp-eq: " << dict["fxp-eq"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "fxp-eq: " << dict["fxp-eq"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  sf64Matrix<D8> fbmul_res(plainA.rows(), plainA.cols());
  start = clock();
  for(int k=0; k<repeats; k++){
      cipher_mul_seq(role, sharedA, lt_res, fbmul_res, eval, enc, runtime);
  }
  end = clock();
  dict["fxp-abmul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "fxp-abmul: " << dict["fxp-abmul"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "fxp-abmul: " << dict["fxp-abmul"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // 4. sint multiplication.
  si64Matrix imul_res(iplainA.rows(), iplainB.cols());
  start = clock();
  for(int k=0; k<repeats; k++)
      eval.asyncMul(runtime, isharedA, isharedB, imul_res).get();
  end = clock();
  dict["int-mul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "int-mul: " << dict["int-mul"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "int-mul: " << dict["int-mul"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // 5. sint addition
  start = clock();
  for(int k=0; k<repeats; k++)
      imul_res = isharedA + isharedB;
  end = clock();
  dict["int-add"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "int-add: " << dict["int-add"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "int-add: " << dict["int-add"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // sint gt
  start = clock();
  for(int k=0; k<repeats; k++)
      cipher_gt(role, isharedB, isharedA, lt_res, eval, runtime);
  end = clock();
  dict["int-gt"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "int-gt: " << dict["int-gt"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "int-gt: " << dict["int-gt"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // sint eq
  start = clock();
  for(int k=0; k<repeats; k++)
      cipher_eq(role, isharedB, isharedA, eq_res, eval, runtime);
  end = clock();
  dict["int-eq"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "int-eq: " << dict["int-eq"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "int-eq: " << dict["int-eq"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  // multiplications
  // 5. sb & si multiplication.
  si64Matrix ibmul_res(iplainA.rows(), iplainA.cols());
  start = clock();
  for(int k=0; k<repeats; k++)
      eval.asyncMul(runtime, isharedA, lt_res, ibmul_res).get();
  end = clock();
  dict["int-abmul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
  for(int r=0; r<task_num; r++){
    if(rank == r){
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "rank: " << rank << std::endl;
      ofs << "int-abmul: " << dict["int-abmul"] << std::endl;  
      ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // if(rank == 0){
  //   // cout << logging_file << endl;
  //   std::ofstream ofs(logging_file, std::ios_base::app);
  //   ofs << "int-abmul: " << dict["int-abmul"] << std::endl;  
  //   ofs.close();
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

    // if(rank == 0){
    //   // cout << logging_file << endl;
    //   std::ofstream ofs(logging_file, std::ios_base::app);
    //   ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime data prep: " << std::setprecision(5) << time_data_prepare << std::endl;

    //   std::map<std::string, double>::iterator iter;
    //   iter = dict.begin();
    //   while(iter != dict.end()){
    //     ofs << iter->first << " " << iter->second << std::endl;
    //     iter ++;
    //   }    
    //   ofs.close();
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

  // }

  return 0;
}