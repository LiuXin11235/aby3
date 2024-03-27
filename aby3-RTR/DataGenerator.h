#include "BuildingBlocks.h"

const long ENCRYPTION_LIMIT = 268435456;

std::vector<aby3::si64> generate_vector_si64(size_t len, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime);
std::vector<int> generate_vector_int(size_t len);