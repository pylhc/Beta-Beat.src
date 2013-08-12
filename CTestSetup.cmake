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

# Find Python..
if(EXISTS "/usr/bin/python2.6") # For our official testing servers, should use this
   set(PYTHON "/usr/bin/python2.6")
else() # For other machines, hope ctest is clever enough..
   find_program(PYTHON python python2)
endif()

