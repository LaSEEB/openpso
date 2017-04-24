cmake_minimum_required(VERSION 2.8.3)
project(SSPSO)

# Compiler options for GCC/Clang
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra -Wall -std=c99")

# Avoid including MinGW dll dependency
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
endif()

add_definitions("-DMT_GENERATE_CODE_IN_HEADER=0")

# Additional include directories
add_subdirectory(libs)

include_directories(libs/mtwist libs/iniparser)

# Create executable
add_executable(${PROJECT_NAME} main.c aloca.c functions.c)
target_link_libraries(${PROJECT_NAME} m mtwist iniparser)
