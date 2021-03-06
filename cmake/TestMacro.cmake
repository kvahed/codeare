macro (MP_TESTS TEST_NAME TEST_CALL)
  if (${WINDOWS})
    set(bat_file ${CMAKE_CURRENT_BINARY_DIR}/test_${TEST_NAME}.bat)
    file(WRITE ${bat_file} "@echo off\n")
    file(APPEND ${bat_file} "SET ${LD_ENV}\n")
    string (REPLACE ";" " " TEST_CALL_STR "${TEST_CALL}")
    file(APPEND ${bat_file} "${TEST_CALL_STR}\n")
    add_test (${TEST_NAME} ${bat_file})
  else ()
    add_test (${TEST_NAME} ${TEST_CALL})
    set_property(TEST ${TEST_NAME} PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${LD_ENV}")
  endif()
endmacro()
