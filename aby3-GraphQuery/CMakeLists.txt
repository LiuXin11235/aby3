
file(GLOB_RECURSE SRCS *.cpp)
add_library(aby3-GraphQuery STATIC ${SRCS})

target_include_directories(aby3-GraphQuery PUBLIC ${CMAKE_SOURCE_DIR}) 
target_link_libraries(aby3-GraphQuery aby3)
