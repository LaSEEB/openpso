# Sources required for building library
set(SRC dictionary.c iniparser.c)

# Add library
add_library(iniparser STATIC ${SRC})
