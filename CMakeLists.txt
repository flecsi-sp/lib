#------------------------------------------------------------------------------#
# Copyright (c) 2023, Triad National Security, LLC.
# All rights reserved.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.23)

#------------------------------------------------------------------------------#
# Module Path
#------------------------------------------------------------------------------#

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#------------------------------------------------------------------------------#
# CMake policies.
#------------------------------------------------------------------------------#

cmake_policy(SET CMP0144 NEW) # ``PACKAGE_NAME''_ROOT
cmake_policy(SET CMP0167 NEW) # FindBoost

#------------------------------------------------------------------------------#
# Project.
#------------------------------------------------------------------------------#

project(FleCSI-SP LANGUAGES NONE)

#------------------------------------------------------------------------------#
# FleCSI-SP Version
#------------------------------------------------------------------------------#

include(version)

#------------------------------------------------------------------------------#
# FleCSI Library
#------------------------------------------------------------------------------#
option(ENABLE_LIBRARY "Enable build of FleCSI-SP library" ON)
mark_as_advanced(ENABLE_LIBRARY)

if(ENABLE_LIBRARY)
  #----------------------------------------------------------------------------#
  # Build shared libraries by default
  #----------------------------------------------------------------------------#
  option(BUILD_SHARED_LIBS "Build shared libs" ON)
  mark_as_advanced(BUILD_SHARED_LIBS)

  enable_language(CXX)
  enable_language(C)
  set(CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_EXTENSIONS OFF)

  set(default_build_type "Debug")
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()

  #------------------------------------------------------------------------------#
  # Find the FleCSI package. This must be done early so that the FleCSI cmake
  # packages will be available to the project, e.g., 'include(FleCSI/format)'.
  #
  # CMake documentation for FleCSI client applications is included in the
  # FleCSI userguide.
  #------------------------------------------------------------------------------#

  if(NOT FleCSI_FOUND)
    find_package(FleCSI 2 REQUIRED)
  endif()

  #----------------------------------------------------------------------------#
  # Define options and includes.
  #----------------------------------------------------------------------------#

  option(ENABLE_EXODUSII "Build with ExodusII Support" ON)

  include(GNUInstallDirs)
  include(FleCSI-SP/exodusii)

  #----------------------------------------------------------------------------#
  # Add library.
  #----------------------------------------------------------------------------#

  add_library(FleCSI-SP)
  add_library(FleCSI-SP::FleCSI-SP ALIAS FleCSI-SP)

  target_include_directories(FleCSI-SP
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_link_libraries(FleCSI-SP PUBLIC FleCSI::FleCSI)

  #----------------------------------------------------------------------------#
  # Add options.
  #----------------------------------------------------------------------------#

  if(ENABLE_EXODUSII)
    flecsi_sp_enable_exodusii(FleCSI-SP)
  endif()

  #----------------------------------------------------------------------------#
  # Capture enable settings.
  #----------------------------------------------------------------------------#

  get_cmake_property(_variableNames VARIABLES)
  string (REGEX MATCHALL "(^|;)ENABLE_[A-Za-z0-9_]*"
    _matchedVars "${_variableNames}")

  foreach(_variableName ${_matchedVars})
    set(FLECSI_${_variableName} ${${_variableName}})
    set(FleCSI_${_variableName} ${${_variableName}})
  endforeach()

  #----------------------------------------------------------------------------#
  # Configure header.
  #----------------------------------------------------------------------------#

  configure_file(${PROJECT_SOURCE_DIR}/flsp/config.hh.in
    ${CMAKE_BINARY_DIR}/flsp/config.hh @ONLY )
    
  install(
    FILES ${CMAKE_BINARY_DIR}/flsp/config.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/flsp
  )

  #------------------------------------------------------------------------------#
  # Extract all project options so they can be exported to the
  # FleCSISPConfig.cmake file.
  #------------------------------------------------------------------------------#

  get_cmake_property(_variableNames VARIABLES)
  string (REGEX MATCHALL "(^|;)FLECSI_[A-Za-z0-9_]*"
    _matchedVars "${_variableNames}")

  foreach(_variableName ${_matchedVars})
    string(REGEX REPLACE "FLSP" "FleCSI-SP" _tmpvariableName ${_variableName})
    set(FLECSI_CONFIG_CODE
      "${FLECSI_CONFIG_CODE}\nset(${_tmpvariableName} \"${${_variableName}}\")")
  endforeach()

  string(REGEX MATCHALL "(^|;)FLOG_[A-Za-z0-9_]*"
    _matchedVars "${_variableNames}")

  foreach(_variableName ${_matchedVars})
    set(FLOG_CONFIG_CODE
      "${FLOG_CONFIG_CODE}\nset(${_variableName} \"${${_variableName}}\")")
  endforeach()

  configure_file(${PROJECT_SOURCE_DIR}/cmake/FleCSI-SP/FleCSISPConfig.cmake.in
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSISPConfig.cmake @ONLY)
  
  install(
    FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSISPConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FleCSI-SP
  )

  #----------------------------------------------------------------------------#
  # Unit Tests
  #----------------------------------------------------------------------------#

  include(FleCSI/unit)

  option(ENABLE_UNIT_TESTS "Enable unit testing" OFF)

  if(ENABLE_UNIT_TESTS)
    set(FLSP_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
    flecsi_enable_testing()
  endif()

  function(flsp_add_unit name)
    flecsi_add_test(${name} ${ARGN})
    if(TARGET ${name})
      target_link_libraries(${name} PUBLIC FleCSI-SP::FleCSI-SP)
    endif()
  endfunction()

  #----------------------------------------------------------------------------#
  # FleCSI-SP library sources.
  #----------------------------------------------------------------------------#

  add_subdirectory(flsp)

  #----------------------------------------------------------------------------#
  # Version
  #----------------------------------------------------------------------------#
  set_target_properties(FleCSI-SP
    PROPERTIES
      VERSION ${FleCSI_SP_VERSION}
  )

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/FleCSSPIConfigVersion.cmake
    VERSION ${FleCSI_SP_VERSION}
    COMPATIBILITY AnyNewerVersion
  )
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FleCSISPConfigVersion.cmake
    DESTINATION
      ${CMAKE_INSTALL_LIBDIR}/cmake/FleCSI-SP
  )
  message("LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
endif()
