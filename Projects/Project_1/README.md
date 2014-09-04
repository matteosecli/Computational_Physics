Project 1
=========
*Solving a simple PDE*

Building:
`g++ main.cpp -o Project_1 -O3 -larmadillo -llapacke -llapack -lgfortran`

Usage:
- `$ ./Project_1 N` solves the PDE using N points with different algorithms and doesn't write output to file (useful for benchmarking purposes).
- `$ ./Project_1 N onealg 1` uses the fastest algorithm and writes output.
- `$ ./Project_1 N onealg 0` uses the fastest algorithm and doesn't write output (useful for testing speed on large N).
