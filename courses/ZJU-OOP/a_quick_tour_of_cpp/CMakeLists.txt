cmake_minimum_required(VERSION 3.10)
project(Sort)
add_executable(Sort test.cpp)

add_custom_target(run COMMAND $<TARGET_FILE:Sort>
		DEPENDS Sort
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_custom_target(myclean COMMAND ${CMAKE_COMMAND} -P
"${CMAKE_BINARY_DIR}/clean.cmake"
		COMMENT "Cleaning build directory")
configure_file("${CMAKE_BINARY_DIR}/clean.cmake"
"${CMAKE_BINARY_DIR}/clean.cmake")
