
# Set binary and source directory..
if("" STREQUAL "${CTEST_SOURCE_DIRECTORY}")
   # Define source directory as current source directory
   file(GLOB _CTEST_SOURCE_DIR CTestTestfile.cmake)
   string(REGEX REPLACE "CTestTestfile.cmake" "" CTEST_SOURCE_DIR ${_CTEST_SOURCE_DIR})
else()
   # Define source directory from run script:
   set(CTEST_SOURCE_DIR ${CTEST_SOURCE_DIRECTORY})
endif()

if(NOT "" STREQUAL "${CTEST_BINARY_DIRECTORY}")
   # Define binary directory from run script:
   set(CTEST_BINARY_DIR ${CTEST_BINARY_DIRECTORY})
endif()
if(NOT DEFINED CTEST_BINARY_DIR)
   # Define binary directory as source directory (not through script)
   set(CTEST_BINARY_DIR ${CTEST_SOURCE_DIR})
endif()

# Define numdiff binary path..
if(WIN32)
   set(NUMDIFF "${CTEST_SOURCE_DIR}/binaries/windows/numdiff-win32.exe")
else() # assume linux, could also be osx..
   set(NUMDIFF "${CTEST_SOURCE_DIR}/binaries/linux/numdiff-linux64")
endif()

# Find Mad-X:
find_program(MADX NAMES madx madx_dev
   HINTS /afs/cern.ch/group/si/slap/bin/
   )

subdirs(GetLLM/test)
subdirs(tests)
