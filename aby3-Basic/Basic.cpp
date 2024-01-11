#include "Basics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>

#include <aby3/sh3/Sh3BinaryEvaluator.h>
#include <aby3/Circuit/CircuitLibrary.h>

#include "../aby3-RTR/debug.h"

using namespace oc;
using namespace aby3;

void bool_cipher_lt(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_lt(bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_eq(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_eq(bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, 1);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_or(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){

    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseOr(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}


void bool_cipher_add(int pIdx, sbMatrix &sharedA, sbMatrix &sharedB, sbMatrix &res, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_add(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}

void bool_cipher_and(int pIdx, aby3::sbMatrix &sharedA, aby3::sbMatrix &sharedB, aby3::sbMatrix &res, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime){
    Sh3BinaryEvaluator binEng;
    CircuitLibrary lib;

    int bitSize = sharedA.bitCount();
    int i64Size = sharedA.i64Size();

    auto cir = lib.int_int_bitwiseAnd(bitSize, bitSize, bitSize);

    binEng.setCir(cir, i64Size, eval.mShareGen);
    binEng.setInput(0, sharedA);
    binEng.setInput(1, sharedB);
    
    auto dep = binEng.asyncEvaluate(runtime).then([&](Sh3Task self) {
        res.resize(i64Size, bitSize);
        binEng.getOutput(0, res);
    });
    dep.get();
}