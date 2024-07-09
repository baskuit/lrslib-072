cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

find_library(GMP_LIBRARY gmp)

set(LIB_SOURCES
    lib.cc
    lrslib.cc
)

# libraries

add_library(lrslib lrsmp.cc ${LIB_SOURCES})
target_compile_definitions(lrslib PRIVATE MP LRS_QUIET)
target_compile_options(lrslib PRIVATE -Wno-write-strings -Wall)

add_library(lrslib64 lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib64 PRIVATE LRSLONG SAFE)
target_compile_options(lrslib64 PRIVATE -Wno-write-strings -Wall)

add_library(lrslib128 lrslong.cc ${LIB_SOURCES})
target_compile_definitions(lrslib128 PRIVATE LRSLONG SAFE B128)
target_compile_options(lrslib128 PRIVATE -Wno-write-strings -Wall)

add_library(lrslibgmp lrsgmp.cc ${LIB_SOURCES})
target_link_libraries(lrslibgmp ${GMP_LIBRARY})
target_compile_definitions(lrslibgmp PRIVATE GMP)
target_compile_options(lrslibgmp PRIVATE -Wno-write-strings -Wall)

# targets

add_executable(main main.cc)
target_link_libraries(main PRIVATE lrslib)

add_executable(main64 main.cc)
target_link_libraries(main64 PRIVATE lrslib64)

add_executable(main128 main.cc)
target_link_libraries(main128 PRIVATE lrslib128)

add_executable(maingmp main.cc)
target_link_libraries(maingmp PRIVATE lrslibgmp)