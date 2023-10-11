#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <fstream>

#ifndef _DEBUG_H_
#define _DEBUG_H_

static std::string debugFile = "/root/aby3/debug.txt";
static std::string debugFolder= "/root/aby3/";

extern void debug_mpi(int rank, int pIdx, std::string info);

extern void debug_info(std::string info);

extern void write_log(std::string log_file, std::string info);

extern void debug_output_vector(std::vector<aby3::si64>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::si64Matrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc);

extern void debug_output_matrix(aby3::sbMatrix& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc, int pIdx, aby3::Sh3Evaluator& eval);

template <aby3::Decimal D>
extern void debug_output_matrix(aby3::sf64Matrix<D>& problem_mat, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){
    aby3::u64 length = problem_mat.rows();    
    aby3::f64Matrix<D> plaininfo(length, 1);
    enc.revealAll(runtime, problem_mat, plaininfo).get();
    // runtime.runUntilTaskCompletes(runtime);

    std::ofstream ofs(debugFile, std::ios_base::app);
    ofs << "length: " << length << std::endl;
    for(int i=0; i<length; i++) ofs << plaininfo(i, 0) << " ";
    ofs << std::endl;
    ofs.close();
}

template <aby3::Decimal D>
extern void debug_output_vector(std::vector<aby3::sf64<D>>& problem_vec, aby3::Sh3Runtime& runtime, aby3::Sh3Encryptor &enc){
    aby3::u64 length = problem_vec.size();
    aby3::sf64Matrix<D> problem_mat; problem_mat.resize(length, 1);
    for(int i=0; i<length; i++) problem_mat(i, 0, problem_vec[i]);
    return debug_output_matrix<D>(problem_mat, runtime, enc); 
}

#endif