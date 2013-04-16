# This is a script for testing the source and submitting
# your results to a common server (cdash).
# You can see the server at:
# http://cern.ch/yngve/cdash/index.php?project=Beta-Beat

# To run, change the source and binary directory to a temporary path
# you have write access to, then run the script with the command:
#   ctest -S run_tests.cmake -VV

## -- SRC Dir
set(CTEST_SOURCE_DIRECTORY "")

## -- BIN Dir
set(CTEST_BINARY_DIRECTORY "")

## -- Dashboard type, possible values are 
##     'Experimental', 'Nightly' or 'Continuous'
set(DASHBOARD Experimental)


## -- Check that the user did what was told:
if(NOT CTEST_SOURCE_DIRECTORY)
    message(FATAL_ERROR "You did not set a source directory")
endif()
if("${CTEST_SOURCE_DIRECTORY}" MATCHES "lintrack")
    message(FATAL_ERROR "You are not allowed to run the tests from lintrack")
endif()

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)

set(CTEST_SITE                          "$ENV{USER}@${HOSTNAME}")

## -- Set site / build name
## --------------------------

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(cpu    -m)

set(CTEST_BUILD_NAME                    "${osname}-${cpu}")




# -----------------------------------------------------------  
# -- commands
# -----------------------------------------------------------  


## -- Update Command
set(CTEST_UPDATE_COMMAND "git")

## -- Checkout command
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
	set(CTEST_CHECKOUT_COMMAND     "${CTEST_UPDATE_COMMAND} clone git://github.com/pylhc/Beta-Beat.src.git ${CTEST_SOURCE_DIRECTORY}")
endif()


ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${DASHBOARD})

## -- Folders containing tests (relative to source directory)
foreach(test_folder . tests GetLLM/test)
   # copy test files to binary folder..
   configure_file("${CTEST_SOURCE_DIRECTORY}/${test_folder}/CTestTestfile.cmake"
                  "${CTEST_BINARY_DIRECTORY}/${test_folder}/CTestTestfile.cmake")
endforeach()

ctest_update()
ctest_test()
ctest_submit()
