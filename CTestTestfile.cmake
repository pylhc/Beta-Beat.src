
# Define source/binary directory from run script:
set(CTEST_SOURCE_DIR ${CTEST_SOURCE_DIRECTORY})
set(CTEST_BINARY_DIR ${CTEST_BINARY_DIRECTORY})

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
