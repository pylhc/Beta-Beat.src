add_test(gen_full_resp_test_output ${PYTHON} ${CTEST_SOURCE_DIR}/MODEL/LHCB/fullresponse/test/test_fullresponse_parallel.py)
set_tests_properties(drive_test_output PROPERTIES TIMEOUT 600)
