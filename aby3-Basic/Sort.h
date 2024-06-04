#include "Basics.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"

int bc_sort_different(std::vector<aby3::sbMatrix> &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size);

int quick_sort_different(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int quick_sort(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int odd_even_merge(aby3::sbMatrix& data1, aby3::sbMatrix& data2, aby3::sbMatrix& res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int odd_even_multi_merge(std::vector<aby3::sbMatrix> &data, aby3::sbMatrix& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);