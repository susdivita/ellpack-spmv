# Project structure

The codebase is split into folders pertaining to:

- [benchmark](./benchmark): defines benchmark suite
- [inputs](./inputs): inputs for the matrices
- [outputs](./outputs): outputs for results (may feed into benchmark)
- [matrix_formats](./matrix_formats): implementations of formats
- [test](./test): testing directory where all tests go
- [build](.build): build results - tipically only used when running commands.

# Utility

- [Google Test](https://google.github.io/googletest/quickstart-cmake.html) - to learn more about how the testing framework works and how to set it up (should work out of the box with the command rules in [CMakeList.txt](./CMakeLists.txt))
- [C++ TestMate](matepek.vscode-catch2-test-adapter) - if you use VSCode, this helps with running tests from the left buttons (super nice interface)

- [launch.json](./VSCode/launch.json) - you must add this when running debugger and tests so that it recognizes the file structures and allows for seemless debugging.

# Building and testing

To build, run from top-level:
`mkdir build`

`cd build && cmake .. && cd ..`

`cmake --build build`

To run tests:

`cd build && ctest && cd ..`

# Running tests

If `TasteMate` got installed, it should recognize the CMake file and automagically recompile and run tests nicely.

# Formats

## [Basic format](./matrix_formats/BasicMatrixFormat.cpp)
## [CSR2](./matrix_formats/CSR2.cpp)
## [BSR](./matrix_formats/BSR.cpp)

# MPI Implementation

Place your implementation for a distributed computation in: [./matrix_format/MPI/](./matrix_formats/MPI/).

Create a new class for each format, which inherits from [./matrix_formats/MPI/AbstractFormatDistributed.h)](./matrix_formats/MPI/AbstractFormatDistributed.h) and implements:
- `InitData()`: initialize the data in the ROOT process
- `ScatterData()`: make sure each process has their window of data
- `PerformMultiplicationDistributed()`: perform MVM

The three methods should be used exactly in this order. At the end, the ROOT process should have the answer vector stored in `ans` (see base class). You can retrieve a C++-vector with `format->GetAnswer();`.

Lastly, update [./matrix_formats/MPI/main.cpp](./matrix_formats/MPI/main.cpp) for your new format.

## Build & run

Build: `make build`

Run: `make run PROCS=2 FORMAT=DDDNaive MTX_IN="../../inputs/Materials\ Problem/arc130.mtx" VEC_IN="../../inputs/Materials\ Problem/test_vec_130_ones_in.txt"`

# Measurements

[LibSciBench](https://spcl.inf.ethz.ch/Research/Performance/LibLSB/)

## Install
- Download the latest version from the link above.
- Go to directory.
- `./configure`
- `make`
- `sudo make install`
- Add installation directory to your `LD_LIBRARY_PATH` variable. (eg. `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib`)
- Compile with `-llsb`.


# Euler setup

## Install LibSciBench
```
wget https://spcl.inf.ethz.ch/Research/Performance/LibLSB/liblsb-0.2.2.tar.gz
tar -xzvf liblsb-0.2.2.tar.gz
cd liblsb-0.2.2
module load intel/19.1.0 openmpi/4.1.4
MPICC=cc MPICXX=CC ./configure --prefix=$(pwd)/build/ --with-mpi --without-papi
make
make install
```

After, add the path to the library to `LD_LIBRARY_PATH`:

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/home/rstiuca/liblsb-0.2.2/build/lib`

## Build
```
module load intel/19.1.0 openmpi/4.1.4
make build_euler
```

## Submit job
See [example of script](./matrix_formats/MPI/euler.sh). See also [Hybrid jobs](https://scicomp.ethz.ch/wiki/Hybrid_jobs).

`sbatch euler.sh`
