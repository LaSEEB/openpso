# Sources required for building library
set(SRC mtwist.c randistrs.c)

add_definitions("-DMT_GENERATE_CODE_IN_HEADER=0")

# Add library
add_library(mtwist SHARED ${SRC})
target_link_libraries(mtwist m)
