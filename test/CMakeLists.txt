# Compile Catch's runner once for all test executables.
add_library(coco_test_main OBJECT main.cpp)
set(coco_tests tests_utils tests_single_encoders_128
    tests_entire_trie_128_big_alphabet tests_entire_trie_128_small_alphabet
    tests_serialize_load tests_boundaries tests_sorted_builder)
foreach(test IN LISTS coco_tests)
    add_executable(${test} $<TARGET_OBJECTS:coco_test_main> ${test}.cpp)
    target_link_libraries(${test} PRIVATE sdsl succinct)
    target_compile_definitions(${test} PRIVATE COCO_DATASET_DIR="${PROJECT_SOURCE_DIR}/dataset")
    if(test MATCHES "big_alphabet|serialize_load|boundaries|sorted_builder")
        target_compile_definitions(${test} PRIVATE BIG_ALPHABET)
    endif()
    add_test(NAME ${test} COMMAND ${test})
endforeach()
