
macro(create_symlink FNAME)
   execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${FNAME} ${CMAKE_CURRENT_BINARY_DIR}/${FNAME})
endmacro()
