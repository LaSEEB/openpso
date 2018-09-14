# Create standalone cec2017_bound_constrained functions library
add_compile_options(
-Wno-unused-result -Wno-maybe-uninitialized -Wno-unused-function)
add_library(funcs_cec2017_bound_constrained SHARED
	cec17_test_func.c cec17_bound_constrained_interface.c)

if (${OPENMP_FOUND} AND ${USE_OPENMP})
	set_target_properties(funcs_cec2017_bound_constrained PROPERTIES
		COMPILE_FLAGS ${OpenMP_C_FLAGS}
		LINK_FLAGS ${OpenMP_C_FLAGS})
endif()
