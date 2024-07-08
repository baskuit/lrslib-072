cmake_minimum_required(VERSION 3.8)
project(lrsnash CXX)

find_library(GMP_LIBRARY gmp)

set(SOURCES
    lrsnash.cc
    lrsnashlib.cc
    lrslib.cc
    lrsdriver.cc
)

set(LIB_SOURCES
    lib.cc
    lrslib.cc
    lrsdriver.cc
)

# old execs

add_executable(lrsnash lrsmp.cc ${SOURCES})
target_compile_definitions(lrsnash PRIVATE MP LRS_QUIET)
target_compile_options(lrsnash PRIVATE -Wno-write-strings -Wall)

add_executable(lrsnash1 lrslong.cc ${SOURCES})
target_compile_definitions(lrsnash1 PRIVATE LRSLONG SAFE)
target_compile_options(lrsnash1 PRIVATE -Wno-write-strings -Wall)

add_executable(lrsnash2 lrslong.cc ${SOURCES})
target_compile_definitions(lrsnash2 PRIVATE LRSLONG SAFE B128)
target_compile_options(lrsnash2 PRIVATE -Wno-write-strings -Wall)

add_executable(lrsnashgmp lrsgmp.cc ${SOURCES})
target_compile_definitions(lrsnashgmp PRIVATE GMP)
target_compile_options(lrsnashgmp PRIVATE -Wno-write-strings -Wall)
target_link_libraries(lrsnashgmp PRIVATE ${GMP_LIBRARY})

# library

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