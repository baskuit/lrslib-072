cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

find_library(GMP_LIBRARY gmp)

set(SOURCES
    lrsnash.cc
    lrsnashlib.cc
    lrslib.cc
    lrsdriver.cc
)

add_executable(lrsnashgmp lrsgmp.cc ${SOURCES})
target_compile_definitions(lrsnashgmp PRIVATE GMP)
target_compile_options(lrsnashgmp PRIVATE -Wno-write-strings -Wall)
target_link_libraries(lrsnashgmp PRIVATE ${GMP_LIBRARY})

add_executable(lrsnash lrsmp.cc ${SOURCES})
target_compile_definitions(lrsnash PRIVATE MP LRS_QUIET)
target_compile_options(lrsnash PRIVATE -Wno-write-strings -Wall)

add_executable(lrsnash1 lrslong.cc ${SOURCES})
target_compile_definitions(lrsnash1 PRIVATE LRSLONG SAFE)
target_compile_options(lrsnash1 PRIVATE -Wno-write-strings -Wall)

add_executable(lrsnash2 lrslong.cc ${SOURCES})
target_compile_definitions(lrsnash2 PRIVATE LRSLONG SAFE B128)
target_compile_options(lrsnash2 PRIVATE -Wno-write-strings -Wall)
