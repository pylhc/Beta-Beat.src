add_test(filecheck_getllm python ${CTEST_SOURCE_DIR}/GetLLM/test/filecheck.py 
				-o 0
				-v ${CTEST_SOURCE_DIR}/GetLLM/test/GetLLM_valid.py
				-m ${CTEST_SOURCE_DIR}/GetLLM/test/../GetLLM.py
				-p ${CTEST_SOURCE_DIR}/GetLLM/test/data
				-s ${CTEST_BINARY_DIR}
				)
set_tests_properties(filecheck_getllm PROPERTIES TIMEOUT 1800)
