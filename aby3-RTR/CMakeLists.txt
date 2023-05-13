
file(GLOB_RECURSE SRCS *.cpp)

add_library(aby3-RTR ${SRCS})

# 判断是否需要使用 MPI 编译  
if(USE_MPI)  
  find_package(MPI REQUIRED)  
  include_directories(${MPI_INCLUDE_PATH})  
  target_link_libraries(aby3-RTR ${MPI_LIBRARIES})  
  set_target_properties(aby3-RTR PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")  
endif()  

target_link_libraries(aby3-RTR aby3)