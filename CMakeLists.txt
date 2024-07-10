cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

set(LIB_SOURCES
    src/lib.cc
    src/lrslib.cc
)

# libraries

add_library(lrslib64 src/lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib64 PRIVATE LRSLONG SAFE)
target_compile_options(lrslib64 PRIVATE -Wno-write-strings -Wall)

add_library(lrslib128 src/lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib128 PRIVATE LRSLONG SAFE B128)
target_compile_options(lrslib128 PRIVATE -Wno-write-strings -Wall)

# targets

add_executable(test64 test.cc)
target_link_libraries(test64 PRIVATE lrslib64)

add_executable(test128 test.cc)
target_link_libraries(test128 PRIVATE lrslib128)

add_executable(benchmark64 benchmark.cc)
target_link_libraries(benchmark64 PRIVATE lrslib64)

add_executable(benchmark128 benchmark.cc)
target_link_libraries(benchmark128 PRIVATE lrslib128)