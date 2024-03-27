#include "Basics.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"

int bc_sort_different(std::vector<aby3::sbMatrix> &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size);

int quick_sort_different(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);