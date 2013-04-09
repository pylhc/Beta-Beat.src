cmake_minimum_required(VERSION 2.6)

project(Beta-Beat C)

enable_testing()
include(CTest)

include(tests/macros.cmake)
add_subdirectory(tests)

add_subdirectory(GetLLM/test)
