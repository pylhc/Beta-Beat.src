add_test(filecheck_getllm python ${CTEST_SOURCE_DIRECTORY}/GetLLM/test/filecheck.py 
				-o 0
				-v ${CTEST_SOURCE_DIRECTORY}/GetLLM/test/GetLLM_valid.py
				-m ${CTEST_SOURCE_DIRECTORY}/GetLLM/test/../GetLLM_dev.py
				-p ${CTEST_SOURCE_DIRECTORY}/GetLLM/test/data
				-s ${CTEST_BINARY_DIRECTORY}
				)
