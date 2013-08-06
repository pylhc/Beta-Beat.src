add_test(getsuper_b1_acdata ${CTEST_SOURCE_DIRECTORY}/tests/getsuper
        ${CTEST_SOURCE_DIRECTORY}
        ${CTEST_SOURCE_DIRECTORY}/tests/models/LHCB1/1_10_1_3
        LHCB1
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@08_21_27_607_0.sdds.new.gz
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@09_40_52_536_0.sdds.new.gz
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam1\@Turn\@2011_08_24\@09_14_23_920_0.sdds.new.gz)

add_test(getsuper_b2_kickdata ${CTEST_SOURCE_DIRECTORY}/tests/getsuper
        ${CTEST_SOURCE_DIRECTORY}
        ${CTEST_SOURCE_DIRECTORY}/tests/models/LHCB2/0500-B2-inj-noac/
        LHCB2
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_27_55_495_0.sdds.new.gz
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_35_57_218_0.sdds.new.gz
        ${CTEST_SOURCE_DIRECTORY}/tests/sdds_files/Beam2\@Turn\@2012_03_15\@09_42_59_188_0.sdds.new.gz)



set_tests_properties(getsuper_b1_acdata getsuper_b2_kickdata PROPERTIES
   TIMEOUT 1200
   )

