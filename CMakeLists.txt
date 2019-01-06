cmake_minimum_required(VERSION 3.8)

project(satellites LANGUAGES CXX CUDA)

enable_language(CUDA)
# set_property(TARGET satellites PROPERTY CUDA_STANDARD 11)

# find_package(CUDA REQUIRED)
# set(CUDA_SEPARABLE_COMPILATION ON)
# set(CUDA_PROPAGATE_HOST_FLAGS OFF)

# # all .cu .cuh files in lib dir
# file(GLOB LIB_CU_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cu)
# file(GLOB LIB_CU_HDRS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cuh)
# # main.cu in crt dir
# file(GLOB MAIN_CU_SRCS main.cu)
# # only for CLion (derps without a cpp file)
# file(GLOB DUMMY_SRCS dummy.cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

file(GLOB srcs "${PROJECT_SOURCE_DIR}/src/*.cpp" 
	           "${PROJECT_SOURCE_DIR}/src/*.h"
	           "${PROJECT_SOURCE_DIR}/src/*.cu"
	           "${PROJECT_SOURCE_DIR}/src/*.cuh")
add_executable(satellites ${srcs})
set_property(TARGET satellites PROPERTY CUDA_STANDARD 11)

# find_package(Threads REQUIRED)
# find_package(OpenMP  REQUIRED)
# target_link_libraries(Threads::Threads)

# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -Wall -Wextra -O2 -march=native -finline-functions -ftree-vectorize")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -Wall -O2 -arch=sm_50")