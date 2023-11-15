cmake_minimum_required(VERSION 3.9)
project(spmmv_dphpc)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(FORMATS_DIRECTORY "./matrix_formats")

set(HEADER_FILES ${FORMATS_DIRECTORY}/AbstractMatrixFormat.h)

enable_testing()

add_library(basicmatrixformat ${FORMATS_DIRECTORY}/BasicMatrixFormat.cpp ${HEADER_FILES})
add_library(csr2 ${FORMATS_DIRECTORY}/CSR2.cpp ${HEADER_FILES})
add_library(csr ${FORMATS_DIRECTORY}/CSR.cpp ${HEADER_FILES})
add_library(bsr ${FORMATS_DIRECTORY}/BSR.cpp ${HEADER_FILES})
add_library(dddnaive ${FORMATS_DIRECTORY}/DDDNaive.cpp ${HEADER_FILES})

add_executable(
  basic_test
  test/all_test.cpp
  ${basicmatrixformat}
  ${csr2}
  ${csr}
  ${bsr}
  ${dddnaive}
  ${HEADER_FILES}
)


target_include_directories(basic_test PRIVATE ${FORMATS_DIRECTORY})


# set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
# set(OPENMP_LIBRARIES "/opt/homebrew/opt/llvm/bin/lib")
# set(OPENMP_INCLUDES "/opt/homebrew/opt/llvm/bin/include")


# Find OpenMP
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
endif()

OPTION (USE_OpenMP "Use OpenMP" ON)

if(USE_OpenMP)
  find_package(OpenMP REQUIRED)
endif(USE_OpenMP)

if (OPENMP_FOUND)
    # include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)
# OpenMP OpenMP
target_link_libraries(basic_test GTest::gtest_main)

# AVX
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSIMDE_ENABLE_OPENMP -fopenmp-simd -O3 -mfma")

include(GoogleTest)
gtest_discover_tests(basic_test)