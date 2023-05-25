#include "PTRTest.h"
#include "PTRFunction.h"

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