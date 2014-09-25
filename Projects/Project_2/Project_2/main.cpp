#include <iostream>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;


// Find the maximum off-diag element
double maxoffdiag ( mat& A, int& k, int& l, int n ) {
    double max = 0.0;
    for ( int i = 0; i < n; i++ ) {
        for ( int j = i + 1; j < n; j++ ) {
            if ( abs(A(i,j)) > max ) {
                max = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
    return max;
}


// This function performs Jacobi rotation
void rotate ( mat& A, mat& R, int k, int l, int n ) {
    // Compute the values of cos and sin
    double s, c;
    if ( A(k,l) != 0.0 ) {
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if ( tau > 0 ) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }

    // Change the elements that have to be changed
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    // Changing the matrix elements with indices k and l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for ( int i = 0; i < n; i++ ) {
        // Changing the remaining elements
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        // Compute the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

/*
Jacobi's method for finding eigenvalues
eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal
of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector
is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n)
n: dimention of matrices
*/
void jacobi_method ( mat& A, mat& R, int n ) {
    // Setting up the eigenvector matrix
    for ( int i = 0; i < n; i++ ) {
        for ( int j = 0; j < n; j++ ) {
            if ( i == j ) {
                R(i,j) = 1.0;
            } else {
                R(i,j) = 0.0;
            }
        }
    }

    // Initializing rotation matrix indexes and constraints
    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = pow(n,3);
    int iterations = 0;
    double max_offdiag = maxoffdiag( A, k, l, n );

    // Keep doing the algorithm while the constraints are respected
    while ( abs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
        max_offdiag = maxoffdiag( A, k, l, n );
        rotate( A, R, k, l, n );
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;
    return;
}

void order_eigenvals(mat& A, mat& R, int n) {
	// Define a workspace double
	double eigenval=0;
	
	// Define a workspace vector
	vec eigenket = zeros<vec>(n-1);
	
	// Define the diagonal elements vector
	vec diag = zeros<vec>(n-1);
    for(int i=0; i<n-1; i++) diag(i) = A(i,i);
	
	// Perform search and mix all the stuff
//	for(int i = 1, i< n-1)
	
	
}





// 'fill_matrix' fills the matrix A in a tridiagonal form.
// Important: matrix A has to be initialized to zero!
void fill_matrix(mat& A, double& rho_min, double& rho_max, int n) {
    // Calculate step lenght
    double h = (rho_max - rho_min) / n;
    double h_square = pow(h,2);

    // Discretize radius and potential
    vec rho(n);
    vec V(n);
    rho(0)=rho_min;
    V(0)=pow(rho_min,2);
    for(int i = 1; i < n; i++){
        rho(i) = rho(i-1) + h;
        V(i) = pow(rho(i),2);
    }

    //Fill the matrix
    A(0,0) = 2.0/h_square + V(1);
    A(0,1) = -1.0/h_square;
    for(int i = 1; i < n-2; i++){
        A(i,i-1) = -1.0/h_square;
        A(i,i) = 2.0/h_square + V(i+1);
        A(i,i+1) = -1.0/h_square;
    }
    A(n-2,n-3) = -1.0/h_square;
    A(n-2,n-2) = 2.0/h_square + V(n-1);
}




int main()
{
    int n = 200;
    mat A = zeros<mat>(n-1,n-1);
    mat R(n-1,n-1);

    double rho_min = 0.0;
    double rho_max = 5.0;

    fill_matrix(A, rho_min, rho_max, n);

//    cout << A << endl;

    jacobi_method(A,R,n-1);
//    cout << A << endl;
//    cout << R << endl;

    ofstream Afile;
    Afile.open("A.txt");
    Afile << A;
    Afile.close();
    ofstream Rfile;
    Rfile.open("R.txt");
    Rfile << R;
    Rfile.close();


    A.reset();
    R.reset();

    return 0;
}

