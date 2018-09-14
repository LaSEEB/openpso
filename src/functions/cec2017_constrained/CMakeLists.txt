# Create standalone cec2017_bound_constrained functions library
add_compile_options(-Wno-unused-result -Wno-unused-parameter)
add_library(funcs_cec2017_constrained SHARED
	cec17_test_COP.c)
#cec15_interface.c)
