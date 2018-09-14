# Create standalone cec2017_bound_constrained functions library
add_compile_options(-Wno-unused-result -Wno-maybe-uninitialized)
add_library(funcs_cec2017_bound_constrained SHARED
	cec17_test_func.c cec17_bound_constrained_interface.c)
