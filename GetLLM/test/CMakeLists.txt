add_test(filecheck_getllm python ${CMAKE_CURRENT_SOURCE_DIR}/filecheck.py 
				-o 0
				-v ${CMAKE_CURRENT_SOURCE_DIR}/GetLLM_valid.py
				-m ${CMAKE_CURRENT_SOURCE_DIR}/../GetLLM.py
				-p ${CMAKE_CURRENT_SOURCE_DIR}/data
				-s ${CMAKE_BINARY_DIR}
				)