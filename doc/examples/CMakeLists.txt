add_info_file(ginac-examples)

add_definitions(-DIN_GINAC)
include_directories(${CMAKE_SOURCE_DIR}/ginac ${CMAKE_BUILD_DIR}/ginac)

add_executable(archive1 archive1.cpp)
target_link_libraries(archive1 ginac ${LIBDL_LIBRARIES})

add_executable(compile1 compile1.cpp)
target_link_libraries(compile1 ginac ${LIBDL_LIBRARIES})

# XXX: compile2 example uses the Cuba library (http://www.feynarts.de/cuba)
# add_executable(compile2 compile2.cpp)
# target_link_libraries(compile2 ginac ${LIBDL_LIBRARIES})

add_executable(compile3 compile3.cpp)
target_link_libraries(compile3 ginac ${LIBDL_LIBRARIES})

add_executable(mystring mystring.cpp)
target_link_libraries(mystring ginac ${LIBDL_LIBRARIES})

add_executable(derivative derivative.cpp)
target_link_libraries(derivative ginac ${LIBDL_LIBRARIES})
