# This is a script for testing the source and submitting
# your results to a common server (cdash).
# You can see the server at:
# http://cern.ch/yngve/cdash/index.php?project=Beta-Beat
# 
# How to:
#  - clone the source from git in a temporary directory:
#       git clone /afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/
#  - set the CTEST_SOURCE_DIRECTORY accordingly below
#  - set a useful test site name (e.g. "your name"."machine type")
#  - set a useful build name (e.g. architecture, os and compiler used)
#  - Run this script with the command:
#          ctest -S cdash_template.cmake

set(CTEST_SITE "myname@maymachine")
set(CTEST_BUILD_NAME "SLC5-64bit")
# Your source should be checked out from git into this directory:
set(CTEST_SOURCE_DIRECTORY "/path/to/bbeat-clone/")
# tests will be ran here:
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/test-folder")

ctest_start(Experimental)

# Do not edit (unless you know cmake/ctest):

set(CTEST_UPDATE_COMMAND "git")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_update()
ctest_configure()
ctest_test()
ctest_submit()
