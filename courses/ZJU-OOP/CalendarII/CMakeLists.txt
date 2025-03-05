cmake_minimum_required(VERSION 3.10)
project(Calendar)
add_executable(Calendar Calendar.cpp)

add_custom_target(run COMMAND $<TARGET_FILE:Calendar>
		DEPENDS Calendar
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR})



