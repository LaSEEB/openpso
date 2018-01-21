cmake_minimum_required(VERSION 3.2)
project(openpso)

# Require C99 without compiler specific extensions
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Compiler options, only GCC and Clang are supported
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}" )

# Avoid including MinGW dll dependency
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
endif()

# Profiling info option
option(PROFILING "Add profiling flags?" OFF)
if (PROFILING)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
	set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg" )
endif()

# Try to use OpenMP
find_package(OpenMP QUIET)
option(USE_OPENMP  "Use OpenMP if available?" ON)

# Definitions required by mtwist library
add_definitions("-DMT_GENERATE_CODE_IN_HEADER=0")

# Directories where to find include files
include_directories(libs/mtwist libs/iniparser libs/zf_log)

# Create standalone functions library
add_library(benchfuncs SHARED src/functions.c src/fprng.c)

# Create standalone PSO library
add_library(${PROJECT_NAME} SHARED src/pso.c
	src/topol/staticgrid2d.c src/topol/staticring1d.c src/topol/staticgraph.c
	libs/mtwist/mtwist.c libs/mtwist/randistrs.c libs/zf_log/zf_log.c)

# Set executable name
set(PSO_EXE_NAME ${PROJECT_NAME}cli)

# Create executable
add_executable(${PSO_EXE_NAME} src/main.c src/pso.c src/functions.c src/fprng.c
	src/topol/staticgrid2d.c src/topol/staticring1d.c src/topol/staticgraph.c
	libs/mtwist/mtwist.c libs/mtwist/randistrs.c libs/iniparser/iniparser.c
	libs/iniparser/dictionary.c libs/zf_log/zf_log.c)
target_link_libraries(${PSO_EXE_NAME} m)

# If OpenMP is found and if we are to use it, set appropriate compiler flags
if (${OPENMP_FOUND} AND ${USE_OPENMP})
	set_target_properties(${PSO_EXE_NAME} PROPERTIES
		COMPILE_FLAGS ${OpenMP_C_FLAGS}
		LINK_FLAGS ${OpenMP_C_FLAGS})
	set_target_properties(${PROJECT_NAME} PROPERTIES
		COMPILE_FLAGS ${OpenMP_C_FLAGS}
		LINK_FLAGS ${OpenMP_C_FLAGS})
	message("Using OpenMP!")
else()
	message("No OpenMP :(")
endif()
