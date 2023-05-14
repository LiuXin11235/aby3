#include "PTRTest.h"
#include "PTRFunction.h"

using namespace oc;
using namespace aby3;
using namespace std;

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