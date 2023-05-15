#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>

static std::string debugFile = "/root/aby3/debug.txt";
static std::string debugFolder= "/root/aby3/";

#ifndef _DEBUG_H_
#define _DEBUG_H_

extern void debug_output_vector(std::vector<aby3::si64>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::si64Matrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::sbMatrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc, int pIdx, aby3::Sh3Evaluator& eval);

extern void debug_mpi(int rank, int pIdx, std::string info);

#endif