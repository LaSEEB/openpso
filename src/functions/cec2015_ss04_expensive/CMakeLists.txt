# Create standalone CEC2015 functions library
add_definitions(-DNO_RECORDING)
add_compile_options(-Wno-unused-result -Wno-maybe-uninitialized)
add_library(funcs_cec2015_ss04_expensive SHARED
	cec15_test_func.c cec15_interface.c)
