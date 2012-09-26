
create_symlink(getsuper)
add_test(getsuper "./getsuper" ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/tests/models/LHCB1/1_10_1_3)

set_tests_properties(getsuper PROPERTIES TIMEOUT 1200)

