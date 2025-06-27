#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find libexoIIv2c
# Find the native EXODUSII headers and libraries.
#
#  EXODUSII_INCLUDE_DIRS - where to find exodusII.h, etc.
#  EXODUSII_LIBRARIES    - List of libraries when using xoIIv2c.
#  EXODUSII_FOUND        - True if exodus found.
#

find_path(EXODUSII_INCLUDE_DIR exodusII.h
  PATH_SUFFIXES include)

#v6.09 calls it libexodus
#debian calls it libexoIIv2
#other distros libexoIIv2c
find_library(EXODUSII_LIBRARY NAMES exodus exoIIv2 exoIIv2c
  PATH_SUFFIXES lib)

set(EXODUSII_LIBRARIES ${EXODUSII_LIBRARY} )
set(EXODUSII_INCLUDE_DIRS ${EXODUSII_INCLUDE_DIR} )

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EXODUSII_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(EXODUSII DEFAULT_MSG EXODUSII_LIBRARY EXODUSII_INCLUDE_DIR )

if(EXODUSII_FOUND AND NETCDF_FOUND)
  set(EXODUSII_LIBRARIES ${EXODUSII_LIBRARY})
  set(EXODUSII_INCLUDE_DIRS ${EXODUSII_INCLUDE_DIR})

  if(NOT TARGET ExodusII::ExodusII)
    add_library(ExodusII::ExodusII UNKNOWN IMPORTED)
    set_target_properties(ExodusII::ExodusII PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${EXODUSII_INCLUDE_DIR}")
    set_property(TARGET ExodusII::ExodusII APPEND PROPERTY
        IMPORTED_LOCATION "${EXODUSII_LIBRARY}")
    target_link_libraries(ExodusII::ExodusII INTERFACE PkgConfig::NETCDF)
  endif()
endif()

mark_as_advanced(EXODUSII_INCLUDE_DIR EXODUSII_LIBRARY)
