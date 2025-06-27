#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find libzoltan
# Find the zoltan headers and libraries.
#
#  EXODUSII_INCLUDE_DIRS - where to find zoltan.h, etc.
#  EXODUSII_LIBRARIES    - List of libraries when using zoltan.
#  EXODUSII_FOUND        - True if exodus found.
#

find_path(ZOLTAN_INCLUDE_DIR zoltan.h
  PATH_SUFFIXES include)

find_library(ZOLTAN_LIBRARY NAMES zoltan
  PATH_SUFFIXES lib)

find_package(ParMETIS REQUIRED QUIET)

set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY} )
set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EXODUSII_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Zoltan DEFAULT_MSG ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR )
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(Zoltan
    REQUIRED_VARS ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR VERSION_VAR ZOLTAN_VERSION_STRING)
else()
  find_package_handle_standard_args(Zoltan
    REQUIRED_VARS ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR)
endif()

if(ZOLTAN_FOUND AND PARMETIS_FOUND)
  set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
  set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR})

  if(NOT TARGET Zoltan::Zoltan)
    add_library(Zoltan::Zoltan UNKNOWN IMPORTED)
    set_target_properties(Zoltan::Zoltan PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${ZOLTAN_INCLUDE_DIR}")
    set_property(TARGET Zoltan::Zoltan APPEND PROPERTY
        IMPORTED_LOCATION "${ZOLTAN_LIBRARY}")
    target_link_libraries(Zoltan::Zoltan INTERFACE ParMetis::ParMetis)
  endif()
endif()

mark_as_advanced(ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY)
