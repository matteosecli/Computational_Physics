#include <iostream>
#include <fstream>
#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;
namespace use {int onealg = 0; int out = 1;}

// 'solvetrid' is a function that solves a linear sistem in 'u' relative to a
// tridiagonal matrix with diagonal elements equal to 'b', subdiagonal elements
// equal to 'a' and superdiagonal elements equal to 'c'. Returns the result
// in the vector 'u', overwriting its elements. The algorithm performs ~8N FLOPS.
void solvetrid(int& N, float& a, float& b, float& c, vec& u, vec& f){
    // Start timing
    clock_t t;
    t = clock();

    // Define the variable 'bet', that is just the denominator of 'gam',
    // and 'gam' itself, that is a workspace vector
    double bet;
    vec gam(N);

    // Start forward substitution
    u[0]=f[0]/(bet=b);
    for(int j = 1; j < N; j++) {
        gam[j]=c/bet;
        bet=b-a*gam[j];
        u[j]=(f[j]-a*u[j-1])/bet;
    }
    // Just one-line backward substitution
    for (int j = (N-2); j >= 0; j--) u[j] -= gam[j+1]*u[j+1];

    // Stop timing and print elapsed time
    t = clock() - t;
    cout << "Elapsed time (solvetrid):\t\t" << ((float)t)/CLOCKS_PER_SEC << "s." << endl;
}

// 'solve_special' is a function that solves a special linear system
// relative to a tridiagonal matrix with b = 2 and a = c = -1. The
// solution has been found analytically, and once the pattern in
// the solution was recognized, it has been coded here. Warning! It
// overwrites 'u' and 'f', so make a copy before calling the function
// if you want to re-use them. The alogrithm performs ~8N FLOPS.
void solve_special(int& N, vec& u, vec& f){
    // Start timing
    clock_t t;
    t = clock();

    for(int j = 1; j < N; j++) f[j] += f[j-1]*j/(j+1);
    u[N-1] = f[N-1]*N/(N+1);
    for(int j = N - 2; j >= 0; j--) u[j] = (f[j] + u[j+1])*(j+1)/(j+2);

    // Stop timing and print elapsed time
    t = clock() - t;
    cout << "Elapsed time (solve_special):\t\t" << ((float)t)/CLOCKS_PER_SEC << "s." << endl;
}

// 'solve_special_idx' is the same as 'solve_special', but it uses
// shifted indexes to reach a (theoretical) minimum of 6N FLOPS.
// However, the number of accesses to the memory is increased, so it
// offers only a small improvement in speed, of the order of 1/1000s
// over 10^8 points.
void solve_special_idx(int& N, vec& u, vec& f){
    // Start timing
    clock_t t;
    t = clock();

    int prev_idx = 1;
    for(int j = 2; j < N + 1; j++) {
        f[j-1] += f[j-2]*prev_idx/j;
        prev_idx = j;
    }
    u[N-1] = f[N-1]*N/(N+1);
    prev_idx = N-1;
    for(int j = N - 1; j > 0; j--) {
        u[j-1] = (f[j-1] + u[j])*j/prev_idx;
        prev_idx = j;
    }

    // Stop timing and print elapsed time
    t = clock() - t;
    cout << "Elapsed time (solve_special_idx):\t" << ((float)t)/CLOCKS_PER_SEC << "s." << endl;
}

// 'split' is a function that discretizes the function 'func',
// storing its values in N points from 0 to 1 in the vector 'f'.
// Grid points are stored in the vector 'x'.
void split(vec& f, vec& x, int& N) {
    // Define points spacing and calculate the grid points
    // and the value of func in those points
    double h = 1.0/(N+1);
    double h_square = pow(h,2);
    for(int i = 0; i < N; i++){
        x[i] = (i+1)*h;
        f[i] = h_square*100*exp(-10*x[i]);
    }
}

// 'relative_error' calculates the relative error with respect to the
// theoretical value 'u_th(x)'.
vec relative_error(vec& u, vec& x, int& N) {
    vec err(N);
    err[0] = 0;
    for(int i = 0; i < N-1; i++){
        err[i] = abs( u[i]/(1-(1-exp(-10))*x[i]-exp(-10*x[i])) - 1 );
    }
    return err;
}

// 'main' takes as first argumt the number of points that the program
// will use during the calculation. Use 'nofile' as second argument
// if you don't want to write output to file (useful for benchmarking
// purposes).
int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);

    // Perform some checks in the optional argument.
    if(argc == 4) {
        use::out = atoi(argv[3]);
        if(strcmp(argv[2], "onealg") == 0 && (strcmp(argv[3], "1") == 0||strcmp(argv[3], "0") == 0)) use::onealg = 1;
        else cout << "Wrong optional argument given. Use 'onealg 1' if you want to use only one algorithm (the fastest) and write to file; use 'onealg 0' if you instead want to write to the output file." << endl;
    }

    // Define the elements of the matrix related to the differential equation
    float a = -1.0;
    float b = 2.0;
    float c = -1.0;

    // Initialize the solution vector 'u' with zeros and the vector 'f'
    // of the function values
    vec u = zeros<vec>(N);
    vec f(N), x(N);
    //for(int i = 0; i < N; i++) f[i] = i+1;

    // Discretize and define workspace vectors
    split(f, x, N);
    vec u_temp(N), f_temp(N), err(N);

    // Compare algorithms only if we want to do benchmarks.
    // This is to save memory if we want just to have grid numbers.
    if(use::onealg == 0) {
        // Solve using 'tridig'
        u_temp = u;
        solvetrid(N, a, b, c, u_temp, f);
        u_temp.reset();

        // Solve using 'solve_special'
        u_temp = u;
        f_temp = f;
        solve_special(N, u_temp, f_temp);
        u_temp.reset();
        f_temp.reset();
    }

    // Solve using 'solve_special_idx'
    solve_special_idx(N, u, f);
    f.reset();

    // Compute relative error
    err = relative_error(u, x, N);
    cout << "Maximum relative error: " << err.max()*100 << "%" << endl;

    // Write the resulting points on the output file
    if(use::onealg == 1 && use::out == 1) {
        // Write the 'x' grid-points to the output file
        ofstream X;
        X.open("X.txt");
        X << x;
        X.close();
        // Write the 'u' grid-points to the output file
        ofstream U;
        U.open("U.txt");
        U << u;
        U.close();
        // Write the error bars to the output file
        ofstream E;
        E.open("E.txt");
        E << err;
        E.close();
    }

    return 0;
}
