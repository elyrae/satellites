cmake_minimum_required(VERSION 3.8)
project(satellites LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
# set(THREADS_PREFER_PTHREAD_FLAG TRUE)

# include_directories(${PROJECT_SOURCE_DIR}/src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
file(GLOB srcs "${PROJECT_SOURCE_DIR}/src/*.cpp" 
	           "${PROJECT_SOURCE_DIR}/src/*.h")
add_executable(satellites ${srcs})
# target_link_libraries( Threads::Threads)

# find_package(Threads REQUIRED)
find_package(OpenMP  REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2 -march=native -finline-functions -ftree-vectorize")