cmake_minimum_required(VERSION 3.8)
project(lrsnash)

find_library(GMP_LIBRARY gmp)

set(GMP_SOURCES
    lrsnashlib.c
    lrslib.c
    lrsgmp.c
    lrsdriver.c
)

set(MP_SOURCES
    lrsnashlib.c
    lrslib.c
    lrsmp.c
    lrsdriver.c
)

set(LONG_SOURCES
    lrsnashlib.c
    lrslib.c
    lrslong.c
    lrsdriver.c
)

add_executable(lrsnashgmp lrsnash.c ${GMP_SOURCES})
target_link_libraries(lrsnashgmp PRIVATE ${GMP_LIBRARY})
target_compile_definitions(lrsnashgmp PRIVATE GMP MA)
set_target_properties(lrsnashgmp PROPERTIES LANG CXX)
set_target_properties(lrsnashgmp PROPERTIES LINKER_LANGUAGE CXX)

add_executable(lrsnash lrsnash.c ${MP_SOURCES})
target_compile_definitions(lrsnash PRIVATE MA MP)
set_target_properties(lrsnash PROPERTIES LANG CXX)
set_target_properties(lrsnash PROPERTIES LINKER_LANGUAGE CXX)

add_executable(lrsnash1 lrsnash.c ${LONG_SOURCES})
target_compile_definitions(lrsnash1 PRIVATE MA LRSLONG SAFE)
set_target_properties(lrsnash1 PROPERTIES LANG CXX)
set_target_properties(lrsnash1 PROPERTIES LINKER_LANGUAGE CXX)

add_executable(lrsnash2 lrsnash.c ${LONG_SOURCES})
target_compile_definitions(lrsnash2 PRIVATE MA LRSLONG SAFE B128)
set_target_properties(lrsnash2 PROPERTIES LANG CXX)
set_target_properties(lrsnash2 PROPERTIES LINKER_LANGUAGE CXX)