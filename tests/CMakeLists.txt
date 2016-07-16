set(BIN_DIR ${CMAKE_BINARY_DIR})
set(TEST_DIR ${CMAKE_BINARY_DIR}/tests)
set(CTS_BIN ${CMAKE_CURRENT_SOURCE_DIR}/cts/bin)
set(CTS_LIB ${CMAKE_CURRENT_SOURCE_DIR}/cts/lib)

add_test(clamr_tests clamr_tests)

add_custom_target(clamr_tests 
                  COMMAND echo "Running clamr tests" &&
                          make clamr_clean_tests &&
                          env BIN_DIR=${BIN_DIR} 
                              CTS_BIN=${CTS_BIN}
                          ${CTS_BIN}/cts -cts clamr.cts clamr.suite &&
                          fgrep SPEEDUP testing/*/*.out
                  DEPENDS ${BIN_DIR}/clamr_cpuonly ${BIN_DIR}/clamr_mpionly
                          ${BIN_DIR}/clamr_mpicheck ${BIN_DIR}/clamr_checkall
                          ${BIN_DIR}/clamr_openmponly ${BIN_DIR}/clamr_mpiopenmponly
                          ${BIN_DIR}/clamr_gpuonly ${BIN_DIR}/clamr_gpucheck
                          ${BIN_DIR}/clamr)

add_custom_target(clamr_clean_tests
                  COMMAND echo "Removing previous tests" && rm -rf ${TEST_DIR}/testing)
