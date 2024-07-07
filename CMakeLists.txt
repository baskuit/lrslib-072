cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

find_library(GMP_LIBRARY gmp)

set(GMP_SOURCES
    lrsnashlib.cc
    lrslib.cc
    lrsgmp.cc
    lrsdriver.cc
)

set(MP_SOURCES
    lrsnashlib.cc
    lrslib.cc
    lrsmp.cc
    lrsdriver.cc
)

set(LONG_SOURCES
    lrsnashlib.cc
    lrslib.cc
    lrslong.cc
    lrsdriver.cc
)

add_executable(lrsnashgmp lrsnash.cc ${GMP_SOURCES})
target_link_libraries(lrsnashgmp PRIVATE ${GMP_LIBRARY})
target_compile_definitions(lrsnashgmp PRIVATE GMP MA)
target_compile_options(lrsnashgmp PRIVATE -Wno-write-strings)

add_executable(lrsnash lrsnash.cc ${MP_SOURCES})
target_compile_definitions(lrsnash PRIVATE MA MP)
target_compile_options(lrsnash PRIVATE -Wno-write-strings)

add_executable(lrsnash1 lrsnash.cc ${LONG_SOURCES})
target_compile_definitions(lrsnash1 PRIVATE MA LRSLONG SAFE)
target_compile_options(lrsnash1 PRIVATE -Wno-write-strings)

add_executable(lrsnash2 lrsnash.cc ${LONG_SOURCES})
target_compile_definitions(lrsnash2 PRIVATE MA LRSLONG SAFE B128)
target_compile_options(lrsnash2 PRIVATE -Wno-write-strings)
