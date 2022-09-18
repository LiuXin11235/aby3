
file(GLOB_RECURSE SRCS *.cpp)

add_library(aby3-RTR ${SRCS})

target_link_libraries(aby3-RTR aby3)