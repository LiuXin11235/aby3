
#include <cryptoTools/Common/CLP.h>
#include <tests_cryptoTools/UnitTests.h>
#include <map>
#include "aby3-Basic/benchmark_basic.h"
#include "eric.h"

using namespace oc;
using namespace aby3;

int main(int argc, char** argv) {
  oc::CLP cmd(argc, argv);

  if(cmd.isSet("BenchmarkORAM")){
    sqrt_oram_benchmark(cmd);
  }

  return 0;
}