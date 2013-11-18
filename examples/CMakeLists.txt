set(EXAMPLES_TESTDATA
   data/cube3072.mesh
   data/cube384.mesh
   data/cube48.mesh
   data/line10.mesh
   data/line23.mesh
   data/mosfet.mesh
   data/nin2d.mesh
   data/square512.mesh
   data/square224.mesh
   data/square128.mesh
   data/square8.mesh
   )
foreach(f IN LISTS EXAMPLES_TESTDATA)
   configure_file(${f} "${CMAKE_CURRENT_BINARY_DIR}/${f}" COPYONLY)
endforeach()

add_subdirectory(tutorial)
