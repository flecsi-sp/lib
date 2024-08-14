function(flsp_headers)
  target_sources(FleCSI-SP PUBLIC FILE_SET public_headers TYPE HEADERS FILES ${ARGN})
endfunction()

function(flsp_sources)
  target_sources(FleCSI-SP PRIVATE ${ARGN})
endfunction()

add_subdirectory(topo)
