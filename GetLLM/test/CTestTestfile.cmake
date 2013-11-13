
file(GLOB tests RELATIVE ${CTEST_SOURCE_DIR}/GetLLM/test/data/ ${CTEST_SOURCE_DIR}/GetLLM/test/data/*)

foreach(test ${tests})
   if(IS_DIRECTORY ${CTEST_SOURCE_DIR}/GetLLM/test/data/${test})
      add_test(getllm_${test} python ${CTEST_SOURCE_DIR}/GetLLM/test/filecheck.py 
          -o 0
          -v ${CTEST_SOURCE_DIR}/GetLLM/test/GetLLM_valid.py
          -m ${CTEST_SOURCE_DIR}/GetLLM/test/../GetLLM.py
          -p ${CTEST_SOURCE_DIR}/GetLLM/test/data
          -s ${CTEST_BINARY_DIR}
          -t ${test}
          )
      set_tests_properties(getllm_${test} PROPERTIES TIMEOUT 900)
   endif()
endforeach()
