# Define numdiff binary path..
if(WIN32)
   set(ENV{NUMDIFF} "${CTEST_SOURCE_DIRECTORY}/binaries/windows/numdiff-win32.exe")
else() # assume linux, could also be osx..
   set(ENV{NUMDIFF} "${CTEST_SOURCE_DIRECTORY}/binaries/linux/numdiff-linux64")
endif()

subdirs(GetLLM/test)
subdirs(tests)
