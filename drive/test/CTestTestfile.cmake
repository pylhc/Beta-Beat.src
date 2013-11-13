add_test(drive_test_output ${PYTHON} ${CTEST_SOURCE_DIR}/drive/test/test_output.py)
set_tests_properties(drive_test_output PROPERTIES TIMEOUT 900)
