# Create standalone CEC2015 functions library
add_compile_options(
-Wno-unused-result -Wno-maybe-uninitialized -Wno-unused-function)
add_library(funcs_cec2015_ss04_expensive SHARED
	cec15_test_func.c cec15_interface.c)

if (${OPENMP_FOUND} AND ${USE_OPENMP})
	set_target_properties(funcs_cec2015_ss04_expensive PROPERTIES
		COMPILE_FLAGS ${OpenMP_C_FLAGS}
		LINK_FLAGS ${OpenMP_C_FLAGS})
endif()
