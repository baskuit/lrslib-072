cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

find_library(GMP_LIBRARY gmp)

set(LIB_SOURCES
    src/lib.cc
    src/lrslib.cc
)

# libraries

add_library(lrslib src/lrsmp.cc ${LIB_SOURCES})
target_compile_definitions(lrslib PRIVATE MP LRS_QUIET)
target_compile_options(lrslib PRIVATE -Wno-write-strings -Wall)

add_library(lrslib64 src/lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib64 PRIVATE LRSLONG SAFE)
target_compile_options(lrslib64 PRIVATE -Wno-write-strings -Wall)

add_library(lrslib128 src/lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib128 PRIVATE LRSLONG SAFE B128)
target_compile_options(lrslib128 PRIVATE -Wno-write-strings -Wall)

add_library(lrslibgmp src/lrsgmp.cc ${LIB_SOURCES})
target_link_libraries(lrslibgmp ${GMP_LIBRARY})
target_compile_definitions(lrslibgmp PRIVATE GMP)
target_compile_options(lrslibgmp PRIVATE -Wno-write-strings -Wall)

# targets

add_executable(test test.cc)
target_link_libraries(test PRIVATE lrslib)

add_executable(test64 test.cc)
target_link_libraries(test64 PRIVATE lrslib64)

add_executable(test128 test.cc)
target_link_libraries(test128 PRIVATE lrslib128)

add_executable(testgmp test.cc)
target_link_libraries(testgmp PRIVATE lrslibgmp)
