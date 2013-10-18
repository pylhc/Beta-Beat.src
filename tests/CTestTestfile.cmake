
macro(getsuper testname accel modelname f1 f2 f3)
   # Test only makes sense if we have found Mad-X and Numdiff:
   if(MADX)
      set(modelpath ${CTEST_SOURCE_DIR}/tests/models/${accel}/${modelname})
      add_test(getsuper_${testname} ${PYTHON}
         ${CTEST_SOURCE_DIR}/GetLLM/getsuper.py
         --madxbin=${MADX}
         --twissfile=${modelpath}/twiss.dat
         --output=./out_${modelname}/ --algorithm=SUSSIX --accel=${accel}
         --beta=${CTEST_SOURCE_DIR}/
         ${f1}
         ${f2}
         ${f3}
         )
      if(NUMDIFF) # Cannot run the numdiff test if the binary has not been found..
      set(NDIFF_PATH ${CTEST_SOURCE_DIR}/tests/${modelname}_out)
         file(GLOB REF_FILES RELATIVE ${NDIFF_PATH}
            ${NDIFF_PATH}/*.out)
         foreach(REF_FILE ${REF_FILES})
            string(REGEX REPLACE ".out" "" FILE_BASE "${REF_FILE}")
            # check if cfg file exists for this output file:
            if(EXISTS ${NDIFF_PATH}/${FILE_BASE}.cfg)
               set(CFG_FILE ${FILE_BASE}.cfg)
            else()
               # select default instead
               set(CFG_FILE default.cfg)
            endif()
            add_test(numdiff_${testname}_${FILE_BASE} ${NUMDIFF}
               -b ${NDIFF_PATH}/${REF_FILE} ./out_${modelname}/${REF_FILE} ${NDIFF_PATH}/${CFG_FILE}
               )
            set_tests_properties(numdiff_${testname}_${REF_FILE}
               PROPERTIES
               FAIL_REGULAR_EXPRESSION "warng"
               DEPENDS getsuper_${testname}
               )
         endforeach()
      endif()
   endif()
endmacro()

# These tests need AFS to run..
if(EXISTS "/afs/cern.ch/")
   getsuper(b1_acdata 
      LHCB1
      1_10_1_3
      ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@08_21_27_607_0.sdds.new.gz
      ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@09_40_52_536_0.sdds.new.gz
      ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@09_14_23_920_0.sdds.new.gz)

   getsuper(b2_kickdata
           LHCB2
           0500-B2-inj-noac
           ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_27_55_495_0.sdds.new.gz
           ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_35_57_218_0.sdds.new.gz
           ${CTEST_SOURCE_DIR}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_42_59_188_0.sdds.new.gz)

   set_tests_properties(getsuper_b1_acdata getsuper_b2_kickdata PROPERTIES
      TIMEOUT 1200
      )
endif()


# test for all PyUnit tests
add_test(all_pyunit_tests python ${CTEST_SOURCE_DIR}/tests/all_tests.py)
set_tests_properties(all_pyunit_tests PROPERTIES TIMEOUT 5400)

