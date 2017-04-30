cmake_minimum_required(VERSION 2.8.3)
project(SSPSO)

# Compiler options for GCC/Clang
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wall -std=c99")

find_package(OpenMP QUIET)

option(USE_OPENMP  "Use OpenMP if available?" ON)

# Avoid including MinGW dll dependency
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
endif()

add_definitions("-DMT_GENERATE_CODE_IN_HEADER=0")

# Additional include directories
add_subdirectory(libs)

include_directories(libs/mtwist libs/iniparser)

# Create functions library
add_library(benchfuncs SHARED functions.c functions_prng.c)

# Create executable
add_executable(${PROJECT_NAME} main.c)
target_link_libraries(${PROJECT_NAME} m mtwist iniparser benchfuncs)

if (${OPENMP_FOUND} AND ${USE_OPENMP})
	set_target_properties(${PROJECT_NAME} PROPERTIES
		COMPILE_FLAGS ${OpenMP_C_FLAGS}
		LINK_FLAGS ${OpenMP_C_FLAGS})
	message("Using OpenMP!")
else()
	message("No OpenMP :(")
endif()
