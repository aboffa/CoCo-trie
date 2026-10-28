cmake_minimum_required(VERSION 3.13)
project(CoCo-trie)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(COCO_NATIVE "Optimize for the build machine's CPU" OFF)
option(COCO_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)
if(COCO_NATIVE)
    add_compile_options(-march=native)
endif()
if(COCO_SANITIZERS)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
endif()

find_package(Boost REQUIRED)
include(cmake/Dependencies.cmake)
include_directories(include "${COCO_DEPS_DIR}/ds2i" "${COCO_DEPS_DIR}/sux")
add_compile_options(-Wall -pedantic)

# Default compaction-depth limits favor small alphabets.
add_executable(coco-trie_bench_small_alphabet benchmark.cpp)
target_link_libraries(coco-trie_bench_small_alphabet sdsl succinct ${Boost_LIBRARIES})

# Shorter compaction-depth limits for larger alphabets.
add_executable(coco-trie_bench_big_alphabet benchmark.cpp)
target_link_libraries(coco-trie_bench_big_alphabet sdsl succinct ${Boost_LIBRARIES})
# The actual alphabet and safe encoding depth are determined from the data.
target_compile_definitions(coco-trie_bench_big_alphabet PRIVATE BIG_ALPHABET)

add_executable(example example.cpp)
target_link_libraries(example sdsl succinct ${Boost_LIBRARIES})
target_compile_definitions(example PRIVATE BIG_ALPHABET)

add_executable(coco-construction-bench scripts/perf/construction/benchmark.cpp)
target_link_libraries(coco-construction-bench PRIVATE sdsl succinct)
target_compile_definitions(coco-construction-bench PRIVATE BIG_ALPHABET)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(test)
endif()
