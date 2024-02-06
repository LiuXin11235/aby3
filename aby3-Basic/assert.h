#include <cstddef>

#ifndef _ASSERT_H_
#define _ASSERT_H_

extern bool checkPowerOfTwo(size_t num) {
    if (num == 0 || (num & (num - 1)) != 0) {
        return false;
    }
    return true;
}

#endif