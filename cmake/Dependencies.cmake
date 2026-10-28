# Upstream CMake files generate files in their source directories. Stage the
# dependencies so both those files and our adaptations stay in the build tree.
set(COCO_DEPS_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps")
function(coco_stage source destination)
    file(GLOB_RECURSE inputs CONFIGURE_DEPENDS "${source}/*")
    list(FILTER inputs EXCLUDE REGEX "/\\.git(/|$)")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${inputs})
    file(COPY "${source}/" DESTINATION "${destination}" PATTERN ".git" EXCLUDE)
endfunction()

file(GLOB ds2i_headers CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/lib/ds2i/*.hpp")
foreach(header IN LISTS ds2i_headers)
    get_filename_component(name "${header}" NAME)
    configure_file("${header}" "${COCO_DEPS_DIR}/ds2i/${name}" COPYONLY)
endforeach()
coco_stage("${PROJECT_SOURCE_DIR}/lib/ds2i/succinct" "${COCO_DEPS_DIR}/ds2i/succinct")
coco_stage("${PROJECT_SOURCE_DIR}/lib/sux/sux" "${COCO_DEPS_DIR}/sux/sux")
coco_stage("${PROJECT_SOURCE_DIR}/lib/sdsl-lite" "${COCO_DEPS_DIR}/sdsl-lite")
# GCC 15 diagnoses these misspelled dependent member names even when the
# affected upstream swap() template is not instantiated.
set(louds_header "${COCO_DEPS_DIR}/sdsl-lite/include/sdsl/louds_tree.hpp")
file(READ "${louds_header}" louds_source)
string(REPLACE "tree.m_select1" "tree.m_bv_select1" louds_source "${louds_source}")
string(REPLACE "tree.m_select0" "tree.m_bv_select0" louds_source "${louds_source}")
file(WRITE "${louds_header}" "${louds_source}")

foreach(header bit_vector broadword intrinsics)
    configure_file("${PROJECT_SOURCE_DIR}/lib/adapted_code/${header}.hpp"
                   "${COCO_DEPS_DIR}/ds2i/succinct/${header}.hpp" COPYONLY)
endforeach()
configure_file("${PROJECT_SOURCE_DIR}/lib/adapted_code/compact_elias_fano.hpp"
               "${COCO_DEPS_DIR}/ds2i/compact_elias_fano.hpp" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/lib/adapted_code/Vector.hpp"
               "${COCO_DEPS_DIR}/sux/sux/util/Vector.hpp" COPYONLY)

set(DS2I_SOURCE_DIR "${COCO_DEPS_DIR}/ds2i")
configure_file("${PROJECT_SOURCE_DIR}/lib/ds2i/ds2i_config.hpp.in"
               "${COCO_DEPS_DIR}/ds2i/ds2i_config.hpp")
set(SUCCINCT_USE_INTRINSICS ON)
configure_file("${PROJECT_SOURCE_DIR}/lib/ds2i/succinct/succinct_config.hpp.in"
               "${COCO_DEPS_DIR}/ds2i/succinct/succinct_config.hpp")
# CoCo only needs the succinct library, not its benchmark executables or Boost
# test libraries. Its public Boost dependency is header-only for this use.
add_library(succinct STATIC
    "${COCO_DEPS_DIR}/ds2i/succinct/rs_bit_vector.cpp"
    "${COCO_DEPS_DIR}/ds2i/succinct/bp_vector.cpp")
target_include_directories(succinct PUBLIC "${COCO_DEPS_DIR}/ds2i" ${Boost_INCLUDE_DIRS})

# Permit the bundled dependencies' older CMake minimums on CMake 4.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
add_subdirectory("${COCO_DEPS_DIR}/sdsl-lite" "${CMAKE_CURRENT_BINARY_DIR}/sdsl" EXCLUDE_FROM_ALL)
