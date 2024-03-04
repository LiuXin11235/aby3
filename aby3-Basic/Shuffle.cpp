#include "Shuffle.h"

using namespace oc;
using namespace aby3;

#define DEBUG_SHUFFLE2
// #define DEBUG_SHUFFLE_COMM

/**
 * The shuffle protocol accoring to
 * https://dl.acm.org/doi/10.1145/3460120.3484560 (advanced version).
 */
int efficient_shuffle(std::vector<sbMatrix>& T, int pIdx,
                      std::vector<sbMatrix>& Tres, Sh3Encryptor& enc,
                      Sh3Evaluator& eval, Sh3Runtime& runtime) {
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.size();
    size_t unit_len = T[0].i64Size();

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);

    // 2 - generate the random masks Z.
    std::vector<i64Matrix> prev_maskZ(len);
    std::vector<i64Matrix> next_maskZ(len);
    for (size_t i = 0; i < len; i++) {
        prev_maskZ[i].resize(unit_len, T[0].bitCount());
        next_maskZ[i].resize(unit_len, T[0].bitCount());
        get_random_mask(pIdx, prev_maskZ[i], prevSeed);
        get_random_mask(pIdx, next_maskZ[i], nextSeed);
    }

    // get the random permutation.
    if (pIdx == 0) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, T[0].bitCount());
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], nextSeed);
            get_random_mask(pIdx, maskA[i], prevSeed);
        }

        // get the sharedX1.
        std::vector<i64Matrix> sharedX1(len);
        for (size_t i = 0; i < len; i++) {
            sharedX1[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX1[i](j, 0) =
                    T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        std::vector<i64Matrix> sharedX2(len);
        for (size_t i = 0; i < len; i++) {
            sharedX2[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedX2);

        // send the sharedX2 to P1.
        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2_matrix(i * unit_len + j, 0) = sharedX2[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedX2_matrix.data(),
                                          sharedX2_matrix.size());

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskA[i](j);
                Tres[i].mShares[0](j, 0) = maskB[i](j);
            }
        }
    }
    if (pIdx == 1) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskB[i], prevSeed);
        }

        // compute the sharedY1 and send to the next party.
        std::vector<i64Matrix> shared_Y1(len);
        for (size_t i = 0; i < len; i++) {
            shared_Y1[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, shared_Y1);

        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedY1_matrix(i * unit_len + j, 0) = shared_Y1[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedY1_matrix.data(),
                                          sharedY1_matrix.size());

        i64Matrix sharedX2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedX2_matrix.data(),
                                 sharedX2_matrix.size());

        // compute the sharedX3.
        std::vector<i64Matrix> sharedX3(len);
        for (size_t i = 0; i < len; i++) {
            sharedX3[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedX3[i](j, 0) =
                    sharedX2_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX3);

        // compute the masked C1.
        // std::vector<i64Matrix> maskedC1(len);
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC1_matrix(i * unit_len + j, 0) =
                    sharedX3[i](j) ^ maskB[i](j);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(maskedC1_matrix.data(),
                                          maskedC1_matrix.size());

        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mNext.recv(maskedC2_matrix.data(),
                                 maskedC2_matrix.size());

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskB[i](j);
                Tres[i].mShares[0](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                           maskedC2_matrix(i * unit_len + j);
            }
        }
    }
    if (pIdx == 2) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskA[i].resize(unit_len, T[0].bitCount());
            get_random_mask(pIdx, maskA[i], nextSeed);
        }

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(sharedY1_matrix.data(),
                                 sharedY1_matrix.size());

        // compute the sharedY2.
        std::vector<i64Matrix> sharedY2(len);
        for (size_t i = 0; i < len; i++) {
            sharedY2[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedY2[i](j, 0) =
                    sharedY1_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedY2);

        std::vector<i64Matrix> sharedY3(len);
        for (size_t i = 0; i < len; i++) {
            sharedY3[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedY3);

        // compute the masked C2.
        i64Matrix maskedC2_matrix(len * unit_len, T[0].bitCount());
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC2_matrix(i * unit_len + j, 0) =
                    sharedY3[i](j) ^ maskA[i](j);
            }
        }
        runtime.mComm.mPrev.asyncSendCopy(maskedC2_matrix.data(),
                                          maskedC2_matrix.size());
        i64Matrix maskedC1_matrix(len * unit_len, T[0].bitCount());
        runtime.mComm.mPrev.recv(maskedC1_matrix.data(),
                                 maskedC1_matrix.size());

        std::vector<i64Matrix> maskC(len);
        for (size_t i = 0; i < len; i++) {
            maskC[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                maskC[i](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                 maskedC2_matrix(i * unit_len + j);
            }
        }

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, T[0].bitCount());
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskC[i](j);
                Tres[i].mShares[0](j, 0) = maskA[i](j);
            }
        }
    }
    return 0;
}

/**
 * The shuffle protocol accoring to
 * https://dl.acm.org/doi/10.1145/3460120.3484560 (advanced version).
 */
int efficient_shuffle_with_random_permutation(
    std::vector<sbMatrix>& T, int pIdx, std::vector<sbMatrix>& Tres,
    std::vector<aby3::si64>& Pi, aby3::Sh3Encryptor& enc,
    aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime) {
    // get the common randomness.
    block prevSeed = enc.mShareGen.mPrevCommon.getSeed();
    block nextSeed = enc.mShareGen.mNextCommon.getSeed();
    size_t len = T.size();
    size_t unit_len = T[0].i64Size();
    size_t bit_count = 1;  // currently, we only support 1-64 bit element representation.

#ifdef DEBUG_SHUFFLE2
    std::ofstream party_fs(PARTY_FILE + std::to_string(pIdx) + ".txt",
                           std::ios::app);
    party_fs << "len = " << len << " unit len = " << unit_len << std::endl;
#endif

    // generate the prev, next - correlated randomness.
    // 1 - generate the permutations and the inverse permutations.
    std::vector<size_t> prev_permutation;
    std::vector<size_t> next_permutation;
    get_permutation(len, prev_permutation, prevSeed);
    get_permutation(len, next_permutation, nextSeed);
    std::vector<size_t> prev_inverse_permutation;
    std::vector<size_t> next_inverse_permutation;
    get_inverse_permutation(prev_permutation, prev_inverse_permutation);
    get_inverse_permutation(next_permutation, next_inverse_permutation);

    // init Pi with the plain shares of R.
    Pi.resize(len);
    for (size_t i = 0; i < len; i++) {
        if (pIdx == 0) {
            Pi[i].mData[0] = ~0;
            Pi[i].mData[0] = i;
        }
        if (pIdx == 1) {
            Pi[i].mData[0] = ~0;
            Pi[i].mData[1] = ~0;
        }
        if (pIdx == 2) {
            Pi[i].mData[0] = i;
            Pi[i].mData[1] = ~0;
        }
    }

    // generate the random masks Z and RZ.
    std::vector<i64Matrix> prev_maskZ(len);
    std::vector<i64Matrix> next_maskZ(len);
    i64Matrix prev_maskRZ(len, bit_count);
    i64Matrix next_maskRZ(len, bit_count);

    for (size_t i = 0; i < len; i++) {
        prev_maskZ[i].resize(unit_len, bit_count);
        next_maskZ[i].resize(unit_len, bit_count);
        get_random_mask(pIdx, prev_maskZ[i], prevSeed);
        get_random_mask(pIdx, next_maskZ[i], nextSeed);
    }
    get_random_mask(pIdx, prev_maskRZ, prevSeed);
    get_random_mask(pIdx, next_maskRZ, nextSeed);

    // get the random permutation.
    if (pIdx == 0) {
        // pre-generate the randomness, masks and the maskRs.
        std::vector<i64Matrix> maskB(len);
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, bit_count);
            maskA[i].resize(unit_len, bit_count);
            get_random_mask(pIdx, maskB[i], nextSeed);
            get_random_mask(pIdx, maskA[i], prevSeed);
        }
        i64Matrix maskRA(len, bit_count);
        get_random_mask(pIdx, maskRA, prevSeed);

        // get the sharedX1.
        std::vector<i64Matrix> sharedX1(len);
        for (size_t i = 0; i < len; i++) {
            sharedX1[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedX1[i](j, 0) =
                    T[i].mShares[0](j) ^ T[i].mShares[1](j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX1);

        // get the sharedX2.
        std::vector<i64Matrix> sharedX2(len);
        for (size_t i = 0; i < len; i++) {
            sharedX2[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2[i](j, 0) = sharedX1[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedX2);

        // send the sharedX2 to P1.
        i64Matrix sharedX2_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedX2_matrix(i * unit_len + j, 0) = sharedX2[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedX2_matrix.data(),
                                          sharedX2_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending X2 to p1, size = " << sharedX2_matrix.size()
                 << std::endl;
#endif

        // the following is for the shared permutation.
        // receive the sharedRY1 from P1.
        i64Matrix sharedRY1_matrix(len, bit_count);

#ifdef DEBUG_SHUFFLE2
        party_fs << "going to receive RY1 from p1, size = " << sharedRY1_matrix.size()
                 << std::endl;
#endif

        runtime.mComm.mNext.recv(sharedRY1_matrix.data(),
                                 sharedRY1_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "received RY1 from p1, size = " << sharedRY1_matrix.size()
                 << std::endl;
#endif

        // compute the sharedRY2.
        std::vector<i64> sharedRY2(len);
        for (size_t i = 0; i < len; i++) {
            sharedRY2[i] = sharedRY1_matrix(i, 0) ^ prev_maskRZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedRY2);

        // compute the sharedRY3.
        std::vector<i64> sharedRY3(len);
        for (size_t i = 0; i < len; i++) {
            sharedRY3[i] = sharedRY2[i] ^ next_maskRZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedRY3);

        i64Matrix maskRB2(len, bit_count);
        for (size_t i = 0; i < len; i++) {
            maskRB2(i, 0) = sharedRY3[i] ^ maskRA(i, 0);
        }
        runtime.mComm.mNext.asyncSendCopy(maskRB2.data(), maskRB2.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending RB2 to p1, size = " << maskRB2.size() << std::endl;
#endif

        i64Matrix maskRB1(len, bit_count);
        runtime.mComm.mNext.recv(maskRB1.data(), maskRB1.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "reseiving RB1 from p1, size = " << maskRB1.size()
                 << std::endl;
#endif

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, BITSIZE);
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskA[i](j);
                Tres[i].mShares[0](j, 0) = maskB[i](j);
            }
        }

        // compute the final shares of Pi.
        for (size_t i = 0; i < len; i++) {
            Pi[i].mData[1] = maskRA(i, 0);
            Pi[i].mData[0] = maskRB1(i, 0) ^ maskRB2(i, 0);
        }
    }
    if (pIdx == 1) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskB(len);
        for (size_t i = 0; i < len; i++) {
            maskB[i].resize(unit_len, bit_count);
            get_random_mask(pIdx, maskB[i], prevSeed);
        }
        i64Matrix maskRC(len, bit_count);
        get_random_mask(pIdx, maskRC, nextSeed);

        // compute the sharedY1 and send to the next party.
        std::vector<i64Matrix> shared_Y1(len);
        for (size_t i = 0; i < len; i++) {
            shared_Y1[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                shared_Y1[i](j, 0) = T[i].mShares[0](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, shared_Y1);

        i64Matrix sharedY1_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                sharedY1_matrix(i * unit_len + j, 0) = shared_Y1[i](j, 0);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(sharedY1_matrix.data(),
                                          sharedY1_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending Y1 to p2, size = " << sharedY1_matrix.size()
                 << std::endl;
#endif

        i64Matrix sharedX2_matrix(len * unit_len, bit_count);

#ifdef DEBUG_SHUFFLE2
        party_fs << "going to receive X2 from p0, size = " << sharedX2_matrix.size()
                 << std::endl;
#endif

        runtime.mComm.mPrev.recv(sharedX2_matrix.data(),
                                 sharedX2_matrix.size());
#ifdef DEBUG_SHUFFLE2
        party_fs << "received X2 from p0, size = " << sharedX2_matrix.size()
                 << std::endl;
#endif
        // compute the sharedX3.
        std::vector<i64Matrix> sharedX3(len);
        for (size_t i = 0; i < len; i++) {
            sharedX3[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedX3[i](j, 0) =
                    sharedX2_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedX3);

        // compute the masked C1.
        // std::vector<i64Matrix> maskedC1(len);
        i64Matrix maskedC1_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC1_matrix(i * unit_len + j, 0) =
                    sharedX3[i](j) ^ maskB[i](j);
            }
        }
        runtime.mComm.mNext.asyncSendCopy(maskedC1_matrix.data(),
                                          maskedC1_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending C1 to p2, size = " << maskedC1_matrix.size()
                 << std::endl;
#endif

        i64Matrix maskedC2_matrix(len * unit_len, bit_count);
        runtime.mComm.mNext.recv(maskedC2_matrix.data(),
                                 maskedC2_matrix.size());


#ifdef DEBUG_SHUFFLE2
        party_fs << "received C2 from p2, size =" << maskedC2_matrix.size()
                 << std::endl;
#endif

        // the following is for shared permutation.
        // compute the sharedRY1
        std::vector<i64> sharedRY1(len);
        for (size_t i = 0; i < len; i++) {
            sharedRY1[i] = Pi[i].mData[1] ^ next_maskRZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedRY1);
        runtime.mComm.mPrev.asyncSendCopy(sharedRY1.data(), sharedRY1.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending RY1 to p0, size = " << sharedRY1.size()
                 << std::endl;
#endif

        // receive the sharedRX2 from the next party.
        i64Matrix sharedRX2_matrix(len, bit_count);
        runtime.mComm.mNext.recv(sharedRX2_matrix.data(),
                                 sharedRX2_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "received RX2 from p2, size = " << sharedRX2_matrix.size()
                 << std::endl;
#endif

        // compute the sharedRX3.
        std::vector<i64> sharedRX3(len);
        for (size_t i = 0; i < len; i++) {
            sharedRX3[i] = sharedRX2_matrix(i, 0) ^ prev_maskRZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedRX3);

        // compute the maskedRB1.
        i64Matrix maskRB1(len, bit_count);
        for (size_t i = 0; i < len; i++) {
            maskRB1(i, 0) = sharedRX3[i] ^ maskRC(i, 0);
        }
        runtime.mComm.mPrev.asyncSendCopy(maskRB1.data(), maskRB1.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending RB1 to p0, size = " << maskRB1.size()
                 << std::endl;
#endif

        i64Matrix maskRB2(len, bit_count);
        runtime.mComm.mPrev.recv(maskRB2.data(), maskRB2.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "received RB2 from p0, size = " << maskRB2.size()
                 << std::endl;
#endif

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, BITSIZE);
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskB[i](j);
                Tres[i].mShares[0](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                           maskedC2_matrix(i * unit_len + j);
            }
        }

        // compute the final shares of Pi.
        for (size_t i = 0; i < len; i++) {
            Pi[i].mData[0] = maskRC(i, 0);
            Pi[i].mData[1] = maskRB1(i, 0) ^ maskRB2(i, 0);
        }
    }
    if (pIdx == 2) {
        // pre-generate the randomness.
        std::vector<i64Matrix> maskA(len);
        for (size_t i = 0; i < len; i++) {
            maskA[i].resize(unit_len, bit_count);
            get_random_mask(pIdx, maskA[i], nextSeed);
        }

        i64Matrix maskRC(len, bit_count);
        i64Matrix maskRA(len, bit_count);
        get_random_mask(pIdx, maskRC, prevSeed);
        get_random_mask(pIdx, maskRA, nextSeed);

        // receive the sharedY1 from the previous party.
        i64Matrix sharedY1_matrix(len * unit_len, bit_count);

#ifdef DEBUG_SHUFFLE2
        party_fs << "going to receive Y1 from p1, size = " << sharedY1_matrix.size()
                 << std::endl;
#endif

        runtime.mComm.mPrev.recv(sharedY1_matrix.data(),
                                 sharedY1_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "received Y1 from p1, size = " << sharedY1_matrix.size()
                 << std::endl;
#endif

        // compute the sharedY2.
        std::vector<i64Matrix> sharedY2(len);
        for (size_t i = 0; i < len; i++) {
            sharedY2[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedY2[i](j, 0) =
                    sharedY1_matrix(i * unit_len + j) ^ next_maskZ[i](j);
            }
        }
        plain_permutate(next_permutation, sharedY2);

        std::vector<i64Matrix> sharedY3(len);
        for (size_t i = 0; i < len; i++) {
            sharedY3[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                sharedY3[i](j, 0) = sharedY2[i](j) ^ prev_maskZ[i](j);
            }
        }
        plain_permutate(prev_permutation, sharedY3);

        // compute the masked C2.
        i64Matrix maskedC2_matrix(len * unit_len, bit_count);
        for (size_t i = 0; i < len; i++) {
            for (size_t j = 0; j < unit_len; j++) {
                maskedC2_matrix(i * unit_len + j, 0) =
                    sharedY3[i](j) ^ maskA[i](j);
            }
        }
        runtime.mComm.mPrev.asyncSendCopy(maskedC2_matrix.data(),
                                          maskedC2_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending C2 to p1, size = " << maskedC2_matrix.size()
                 << std::endl;
#endif

        i64Matrix maskedC1_matrix(len * unit_len, bit_count);
        runtime.mComm.mPrev.recv(maskedC1_matrix.data(),
                                 maskedC1_matrix.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "received C1 from p1, size = " << maskedC1_matrix.size()
                 << std::endl;
#endif

        std::vector<i64Matrix> maskC(len);
        for (size_t i = 0; i < len; i++) {
            maskC[i].resize(unit_len, bit_count);
            for (size_t j = 0; j < unit_len; j++) {
                maskC[i](j, 0) = maskedC1_matrix(i * unit_len + j) ^
                                 maskedC2_matrix(i * unit_len + j);
            }
        }

        // the following is for the shared permutation.

        // compute the sharedRX1.
        std::vector<i64> sharedRX1(len);
        for (size_t i = 0; i < len; i++) {
            sharedRX1[i] = Pi[i].mData[0] ^ Pi[i].mData[1] ^ prev_maskRZ(i, 0);
        }
        plain_permutate(prev_inverse_permutation, sharedRX1);

        // compute the sharedRX2.
        std::vector<i64> sharedRX2(len);
        for (size_t i = 0; i < len; i++) {
            sharedRX2[i] = sharedRX1[i] ^ next_maskRZ(i, 0);
        }
        plain_permutate(next_inverse_permutation, sharedRX2);

        runtime.mComm.mPrev.asyncSendCopy(sharedRX2.data(), sharedRX2.size());

#ifdef DEBUG_SHUFFLE2
        party_fs << "sending RX2 to p1, size = " << sharedRX2.size()
                 << std::endl;
#endif

        // compute the final shares.
        for (size_t i = 0; i < len; i++) {
            Tres[i].resize(unit_len, BITSIZE);
            for (size_t j = 0; j < unit_len; j++) {
                Tres[i].mShares[1](j, 0) = maskC[i](j);
                Tres[i].mShares[0](j, 0) = maskA[i](j);
            }
        }

        // compute the final shares of Pi.
        for (size_t i = 0; i < len; i++) {
            Pi[i].mData[0] = maskRA(i, 0);
            Pi[i].mData[1] = maskRC(i, 0);
        }
    }

    return 0;
}