#include "PTRTest.h"
#include "PTRFunction.h"
#include "BuildingBlocks.h"

using namespace oc;
using namespace aby3;
using namespace std;

// #define DEBUG

int test_cipher_index_ptr(CLP& cmd, int n, int m){

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


int test_cipher_index_ptr_mpi(CLP& cmd, int n, int m, int task_num, int opt_B){

  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout << rank << endl;

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
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(vecIndex.data(), range_index, vecM.data(), res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "time_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }

  // // // #ifdef DEBUG
  // if(rank == 0){
  //   cout << "my rank = " << rank << " | I AM IN HERE!!!" << endl;
  //   // test for correctness:
  //   si64Matrix aby3res(m, 1);
  //   for(int i=0; i<m; i++){
  //     aby3res(i, 0, res[i]);
  //   }

  //   i64Matrix plainTest(n, 1);
  //   for (int i = 0; i < n; i++) plainTest(i, 0) = i;

  //   i64Matrix plain_res;
  //   enc.revealAll(runtime, aby3res, plain_res).get();
  //   if (role == 0) {
  //     for(int i=0; i<m; i++){
  //       if(plain_res(i, 0) != plainTest(n - 1 - i)){
  //         cout << "res != test in index " << i << "\nres[i] = " <<  plain_res(i, 0) << " plainTest = " << plainTest(n - 1 - i) << endl;
  //       }
  //     }
  //     cout << "FINISH!" << endl;
  //   }
  // }
// #endif

}


int test_cipher_select_ptr_mpi(CLP& cmd, int n, int m, int task_num, int opt_B){

  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout << rank << endl;

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
  // evaluate the task.
  mpiPtrTask->circuit_evaluate(range_index, vecIndex.data(), vecM.data(), res.data());
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "time_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }

  // // // #ifdef DEBUG
  // if(rank == 0){
  //   cout << "my rank = " << rank << " | I AM IN HERE!!!" << endl;
  //   // test for correctness:
  //   si64Matrix aby3res(m, 1);
  //   for(int i=0; i<m; i++){
  //     aby3res(i, 0, res[i]);
  //   }

  //   i64Matrix plainTest(n, 1);
  //   for (int i = 0; i < n; i++) plainTest(i, 0) = i;

  //   i64Matrix plain_res;
  //   enc.revealAll(runtime, aby3res, plain_res).get();
  //   if (role == 0) {
  //     for(int i=0; i<m; i++){
  //       if(plain_res(i, 0) != plainTest(n - 1 - i)){
  //         cout << "res != test in index " << i << "\nres[i] = " <<  plain_res(i, 0) << " plainTest = " << plainTest(n - 1 - i) << endl;
  //       }
  //     }
  //     cout << "FINISH!" << endl;
  //   }
  // }
// #endif

}


int test_cipher_rank_ptr_mpi(oc::CLP& cmd, int n, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout << rank << endl;

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
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "time_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }

  // // // #ifdef DEBUG
  // if(rank == 0){
  //   cout << "my rank = " << rank << " | I AM IN HERE!!!" << endl;
  //   // test for correctness:
  //   si64Matrix aby3res(m, 1);
  //   for(int i=0; i<m; i++){
  //     aby3res(i, 0, res[i]);
  //   }

  //   i64Matrix plainTest(n, 1);
  //   for (int i = 0; i < n; i++) plainTest(i, 0) = i;

  //   i64Matrix plain_res;
  //   enc.revealAll(runtime, aby3res, plain_res).get();
  //   if (role == 0) {
  //     for(int i=0; i<m; i++){
  //       if(plain_res(i, 0) != plainTest(n - 1 - i)){
  //         cout << "res != test in index " << i << "\nres[i] = " <<  plain_res(i, 0) << " plainTest = " << plainTest(n - 1 - i) << endl;
  //       }
  //     }
  //     cout << "FINISH!" << endl;
  //   }
  // }
// #endif
}


int test_cipher_sort_ptr_mpi(oc::CLP& cmd, int n, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout << rank << endl;

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
  end = clock();
  double time_task_eval = double((end - start)*1000)/(CLOCKS_PER_SEC);

  if(rank == 0){
    // cout << logging_file << endl;
    std::ofstream ofs(logging_file, std::ios_base::app);
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "time_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }

  // // // #ifdef DEBUG
  // if(rank == 0){
  //   cout << "my rank = " << rank << " | I AM IN HERE!!!" << endl;
  //   // test for correctness:
  //   si64Matrix aby3res(m, 1);
  //   for(int i=0; i<m; i++){
  //     aby3res(i, 0, res[i]);
  //   }

  //   i64Matrix plainTest(n, 1);
  //   for (int i = 0; i < n; i++) plainTest(i, 0) = i;

  //   i64Matrix plain_res;
  //   enc.revealAll(runtime, aby3res, plain_res).get();
  //   if (role == 0) {
  //     for(int i=0; i<m; i++){
  //       if(plain_res(i, 0) != plainTest(n - 1 - i)){
  //         cout << "res != test in index " << i << "\nres[i] = " <<  plain_res(i, 0) << " plainTest = " << plainTest(n - 1 - i) << endl;
  //       }
  //     }
  //     cout << "FINISH!" << endl;
  //   }
  // }
// #endif
}


int test_cipher_search_ptr_mpi(oc::CLP& cmd, int n, int m, int task_num, int opt_B){
  // Get current process rank and size  
	int rank, size;  
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  cout << rank << endl;

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
    ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "time_data_prepare: " << std::setprecision(5) << time_task_prep << "\ntime_task_init: " << std::setprecision(5) << time_task_init << "\ntime_task_evaluate: " << std::setprecision(5) << time_task_eval << "\n" << std::endl;
    ofs.close();
  }

  return 0;
}


int test_vectorization(oc::CLP& cmd, int n_, int task_num){
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

  vector<int> n_list = {10,       20,       42,       88,      183,      379,
            784,     1623,     3359,     6951,    14384,    29763,
          61584,   127427,   263665,   545559,  1128837,  2335721,
        4832930, 10000000};
  vector<int> repeats_list = {500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 100, 100, 100, 50, 50, 50, 50, 50, 10, 10};

  // vector<int> n_list = {
  //           784,     1623,     3359,     6951,    14384,    29763,
  //         61584,   127427,   263665,   545559,  1128837,  2335721,
  //       4832930, 10000000};
  // vector<int> repeats_list = {500, 500, 500, 500, 100, 100, 100, 50, 50, 50, 50, 50, 10, 10};

  // vector<int> n_list = {3359,     6951,    14384,    29763,
  //         61584,   127427,   263665,   545559,  1128837,  2335721,
  //       4832930, 10000000};
  // vector<int> repeats_list = {500, 500, 100, 100, 100, 50, 50, 50, 50, 50, 10, 10};

  for(int p = 0; p<n_list.size(); p++){

    int n = n_list[p];
    int repeats = repeats_list[p];

    std::string logging_file = LOG_FOLDER + "log-config-N=" + std::to_string(n) + "-TASKS=" + std::to_string(task_num) + "-" + std::to_string(rank);

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
    MPI_Barrier(MPI_COMM_WORLD);

    // begin test functions.
    // 1. sfixed multiplication.
    sf64Matrix<D8> mul_res(plainA.rows(), plainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        eval.asyncMul(runtime, sharedA, sharedB, mul_res).get();
    end = clock();
    dict["fxp-mul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // 2. sfixed addition
    sf64Matrix<D8> add_res(plainA.rows(), plainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        add_res = sharedA + sharedB;
    end = clock();
    dict["fxp-add"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // 3. sfixed greater-then
    sbMatrix lt_res;
    start = clock();
    for(int k=0; k<repeats; k++)
        cipher_gt(role, sharedB, sharedA, lt_res, eval, runtime);
    end = clock();
    dict["fxp-gt"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    sbMatrix eq_res;
    start = clock();
    for(int k=0; k<repeats; k++)
        cipher_eq(role, sharedB, sharedA, eq_res, eval, runtime);
    end = clock();
    dict["fxp-eq"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    sf64Matrix<D8> fbmul_res(plainA.rows(), plainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++){
        cipher_mul_seq(role, sharedA, lt_res, fbmul_res, eval, enc, runtime);
    }
    end = clock();
    dict["fxp-abmul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // 4. sint multiplication.
    si64Matrix imul_res(iplainA.rows(), iplainB.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        eval.asyncMul(runtime, isharedA, isharedB, imul_res).get();
    end = clock();
    dict["int-mul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // 5. sint addition
    start = clock();
    for(int k=0; k<repeats; k++)
        imul_res = isharedA + isharedB;
    end = clock();
    dict["int-add"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // sint gt
    start = clock();
    for(int k=0; k<repeats; k++)
        cipher_gt(role, isharedB, isharedA, lt_res, eval, runtime);
    end = clock();
    dict["int-gt"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // sint eq
    start = clock();
    for(int k=0; k<repeats; k++)
        cipher_eq(role, isharedB, isharedA, eq_res, eval, runtime);
    end = clock();
    dict["int-eq"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    // multiplications
    // 5. sb & si multiplication.
    si64Matrix ibmul_res(iplainA.rows(), iplainA.cols());
    start = clock();
    for(int k=0; k<repeats; k++)
        eval.asyncMul(runtime, isharedA, lt_res, ibmul_res).get();
    end = clock();
    dict["int-abmul"] = double((end - start)*1000)/(CLOCKS_PER_SEC * repeats);
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
      // cout << logging_file << endl;
      std::ofstream ofs(logging_file, std::ios_base::app);
      ofs << "time_setup: " << std::setprecision(5) << time_task_setup << "\ntime data prep: " << std::setprecision(5) << time_data_prepare << std::endl;

      std::map<std::string, double>::iterator iter;
      iter = dict.begin();
      while(iter != dict.end()){
        ofs << iter->first << " " << iter->second << std::endl;
        iter ++;
      }    
      ofs.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);

  }

  return 0;
}