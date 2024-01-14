
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3_tests/Test.h"
#include "aby3_tests/aby3_tests.h"
#include "eric.h"


using namespace oc;
using namespace aby3;
// std::vector<std::string> unitTestTag{ "u", "unitTest" };

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);
  // reinit the environment and then finalize the environment.
//   MPI_Init(&argc, &argv);

  // if (cmd.isSet(unitTestTag)){
  // 	auto tests = aby3_tests;
  //   tests.runIf(cmd);
  //   return 0;
  // }

  if (cmd.isSet("Bool")){
    bool_basic_test(cmd);
  }

  if (cmd.isSet("Arith")){
    arith_basic_test(cmd);
  }

  return 0;
}