#include "Basics.h"

static const size_t MAX_SENDING_SIZE = 1<<27;

int large_data_sending(int pIdx, aby3::i64Matrix &sharedA, aby3::Sh3Runtime &runtime, bool toNext){
    size_t len = sharedA.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    if(pIdx == 0) debug_info("In sending - len = " + std::to_string(len) + ", round = " + std::to_string(round) + ", last_len = " + std::to_string(last_len));

    if(round == 1){
        if(toNext) runtime.mComm.mNext.asyncSend(sharedA.data(), sharedA.size());
        else runtime.mComm.mPrev.asyncSend(sharedA.data(), sharedA.size());
        return 0;
    }

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix sending_data = sharedA.block(i * MAX_SENDING_SIZE, 0, sending_len, 1);
        if(toNext) runtime.mComm.mNext.asyncSend(sending_data.data(), sending_data.size());
        else runtime.mComm.mPrev.asyncSend(sending_data.data(), sending_data.size());
    }

    return 0;
}

int large_data_receiving(int pIdx, aby3::i64Matrix &res, aby3::Sh3Runtime &runtime, bool fromPrev){
    size_t len = res.rows();
    size_t round = (size_t)ceil(len / (double)MAX_SENDING_SIZE);
    size_t last_len = len - (round - 1) * MAX_SENDING_SIZE;

    if(pIdx == 0) debug_info("In receving - len = " + std::to_string(len) + ", round = " + std::to_string(round) + ", last_len = " + std::to_string(last_len));

    if(round == 1){
        if(fromPrev) runtime.mComm.mPrev.recv(res.data(), res.size());
        else runtime.mComm.mNext.recv(res.data(), res.size());
        return 0;
    }

    for(size_t i=0; i<round; i++){
        size_t sending_len = (i == round - 1) ? last_len : MAX_SENDING_SIZE;
        aby3::i64Matrix receiving_data(sending_len, 1);
        if(fromPrev) runtime.mComm.mPrev.recv(receiving_data.data(), receiving_data.size());
        else runtime.mComm.mNext.recv(receiving_data.data(), receiving_data.size());
        res.block(i * MAX_SENDING_SIZE, 0, sending_len, 1) = receiving_data;
    }

    return 0;
}