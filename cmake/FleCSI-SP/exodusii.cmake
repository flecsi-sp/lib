macro(flecsi_sp_enable_exodusii target)
  find_package(EXODUSII REQUIRED)
  message(STATUS "Found EXODUSII: ${EXODUSII_INCLUDE_DIR}")

  target_include_directories(${target} SYSTEM INTERFACE ${EXODUSII_INCLUDE_DIR})
  target_link_libraries(${target} INTERFACE ${EXODUSII_LIBRARY})
endmacro()
