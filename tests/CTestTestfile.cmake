add_test(getsuper ${CTEST_SOURCE_DIRECTORY}/tests/getsuper
		  ${CTEST_SOURCE_DIRECTORY}
		  ${CTEST_SOURCE_DIRECTORY}/tests/models/LHCB1/1_10_1_3)

set_tests_properties(getsuper PROPERTIES
   TIMEOUT 1200
   FAIL_REGULAR_EXPRESSION "(Value|Type|ZeroDivision|Name)Error"
   )

