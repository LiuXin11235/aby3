#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include <mpi.h>
#include "aby3-GraphQuery/benchmark.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);

  if(cmd.isSet("privGraph")){
    privGraph_performance_profiling(cmd);
  }
  if(cmd.isSet("adjmat")){
    adj_performance_profiling(cmd);
  }
  if(cmd.isSet("edgelist")){
    list_performance_profiling(cmd);
  }

  return 0;
}