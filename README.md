### Summary

_OpenPSO_ is an efficient, modular and multicore-aware framework for
experimenting with different PSO approaches. The package is implemented in C99,
and transparently parallelized with OpenMP. _OpenPSO_ is composed of three
modules:

1. A PSO algorithm library, with the capability of performing parallel function
evaluations.
2. A library of benchmarking functions.
3. A command-line tool for directly experimenting with the different PSO
algorithms and benchmarking functions.

The library components can be interfaced with other programs and programming
languages, making _OpenPSO_ a flexible and adaptable framework for PSO research.

### Building

_OpenPSO_ has been tested with GCC and Clang on Windows and Linux. OpenMP
support depends on the compiler version, but _OpenPSO_ will work regardless
(albeit slower without OpenMP support).  

_OpenPSO_ uses the [CMake] build system, which is able to generate projects for
different targets, e.g. regular Makefiles, XCode or Eclipse.

A [Dev-C++] project file is included for convenience. This project will only
generate the command-line tool. It will not generate the library components for
interfacing with third-party applications. Unfortunately, due to a bug in
[Dev-C++], the editor crashes when opening the project. To avoid this, disable
"Enable code completion" in Tools=>Editor Options=>Completion before opening
the project.

### Experimenting

The `openpsocli` tool uses the following syntax:

`openpsocli [INPUT_FILE [SEED]]`

If no input file is given, the tool defaults to `input.ini`. The input file
defines the PSO parameters and number of runs to perform. An example input file
is available [here](input.ini).

By default, the number of threads used is the same as the number of available
processors. However, this default can be overridden by setting the
`OMP_NUM_THREADS` environment to a specific value.

### License

[Mozilla Public License 2.0](LICENSE)

<!-- links -->
[CMake]: https://cmake.org/
[Dev-C++]: http://orwelldevcpp.blogspot.pt/
