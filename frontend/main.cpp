
#include <cryptoTools/Common/CLP.h>
#include <map>
#include <mpi.h>
#include "aby3-RTR/Test.h"
#include "eric.h"


using namespace oc;
using namespace aby3;


int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);
  // reinit the environment and then finalize the environment.
  MPI_Init(&argc, &argv);

  if (cmd.isSet("Test")){
    basic_test(cmd);
    return 0;
  }

}