
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3-RTR/PtATests.h"
#include "eric.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);

  MPI_Init(&argc, &argv);

  if(cmd.isSet("cipher_index")) {
	test_cipher_index_pta(cmd);
  }

  MPI_Finalize();

  return 0;
}