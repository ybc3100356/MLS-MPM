# MLS-MPM

This is the repo for our implementation of the paper A Moving Least Squares Material Point Method with Displacement Discontinuity and Two-Way Rigid Body Coupling, ACM Transactions on Graphics (SIGGRAPH 2018).

[paper link](https://dl.acm.org/doi/10.1145/3197517.3201293)

## Build

The code is tested on RTX 2060 and RTX 3070 with real-time performance on Linux.

Our implemenation uses CUDA to make the simulation real-time. Please make sure you have CUDA SDK installed.

Entering any of the projects below, you can build them easily using CMake.

```bash
mkdir build
cd build
cmake ..
make -j8
```

## Run

After building, you should see a binary in the build folder. Run that binary and should see the result.
