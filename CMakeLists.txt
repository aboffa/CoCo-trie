cmake_minimum_required(VERSION 3.13)
project(CoCo-trie)

set(CMAKE_CXX_STANDARD 20)

configure_file(
        lib/ds2i/ds2i_config.hpp.in
        lib/ds2i/ds2i_config.hpp
        ESCAPE_QUOTES)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -pedantic")

find_package(Boost REQUIRED)

include_directories(include lib/ds2i lib/sux)

add_subdirectory(lib/ds2i/succinct)
add_subdirectory(lib/sdsl-lite)

#alphabets at most 26
add_executable(coco-trie_bench_small_alphabet benchmark.cpp)
target_link_libraries(coco-trie_bench_small_alphabet sdsl succinct ${Boost_LIBRARIES})

#alphabets at most 95
add_executable(coco-trie_bench_big_alphabet benchmark.cpp)
target_link_libraries(coco-trie_bench_big_alphabet sdsl succinct ${Boost_LIBRARIES})
# this line sets the maximum alphabet size. By default alphabet size < 26
target_compile_definitions(coco-trie_bench_big_alphabet PRIVATE BIG_ALPHABET)

add_executable(example example.cpp)
target_link_libraries(example sdsl succinct ${Boost_LIBRARIES})
target_compile_definitions(example PRIVATE BIG_ALPHABET)

enable_testing()
add_subdirectory(test)