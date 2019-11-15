# Create standalone basic functions library
add_library(funcs_basic SHARED
	functions.c ../../../libs/mtwist/mtwist.c ../../../libs/mtwist/randistrs.c)
