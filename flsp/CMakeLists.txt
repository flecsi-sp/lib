function(flsp_headers)
  target_sources(FleCSI-SP INTERFACE FILE_SET public_headers TYPE HEADERS FILES ${ARGN})
endfunction()

flsp_headers(
  unstructured.hh
)

add_subdirectory(unstructured)
