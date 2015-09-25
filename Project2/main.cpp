/*
 PROJECT 2:
 Jacobi's method for finding eigenvalues and eigenvectors of the symmetric matrix A.

 A: input matrix (n x n)
 R: empty matrix for eigenvectors (n x n)
 n: dimension of matrices
 */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;



// FUNCTIONS

// Calculates and return the value of the potential for a given value of x.
double potential(double x){
    return x*x; // harmonic oscillator potential
}

// Calculates and return the value of the potential for a given value of x and omega.
double potential_coulomb(double x, double omega){
    return omega*omega*x*x + 1./x; // harmonic oscillator potential with coulomb interaction
}


// Find maximum matrix element
double maxoffdiag (double ** A, int * k, int * l, int n)
{
    double max = 0.0;

    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j ++){
            if (fabs(A[i][j]) > max){
                max = fabs(A[i][j]); // value of max element
                *l = i; // row index of max element
                *k = j; // column index of max element
            }
        }
    }
    return max;
}


// Find values of cos and sin
void rotate (double ** A, double ** R, int k, int l, int n)
{
    double s,c; // s = sin(theta), c = cos(theta)

    if (A[k][l] != 0.0){
        double t, tau; // t = tan(theta), tau = cot(2*theta)
        tau = (A[l][l] - A[k][k]/(2*A[k][l]));

        if (tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = - 1.0/(- tau + sqrt(1.0 + tau*tau));
        }
        c = 1./sqrt(1 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];

    // Changing the matrix elements with indices k and l
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + s*s*a_ll;
    A[k][l] = 0.0;
    A[l][k] = 0.0;

    // Changing the remaining elements
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }
        // Compute new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }
    return;
}


// Compute the eigenvalues ??
//void jacobi_method (double ** A, double ** R, int n){
//    // something something
//    return;
//}




// MAIN PROGRAM
int main()
{
    clock_t start, finish; // declare start and final time
    start = clock();

    int n = 10;
    int k,l;
    double rho_min = 0.0; // minimum value
    double rho_max = 10.0; // maximum value (set by the user)
    double h = (rho_max - rho_min)/( (double) n); // step length

    // Set up eigenvector matrix
    double ** R;
    R = new double*[n];
    for (int i = 0; i < n; i++)
        R[i] = new double[n];

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                R[i][j] = 1.0; // diagonal matrix elements
            }
            else{
                R[i][j] = 0.0; // non-diagonal matrix elements
            }
        }
    }

    cout << R << endl;


    // Set up potential
    double rho_i,V[n];

    for (int i = 0; i < n; i++){
        rho_i = rho_min + i*h;
        V[i] = potential(rho_i);
    }


    // Set up tridiagonal matrix
    double ** A;
    A = new double*[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    double e = -1./(h*h); // non-diagonal constant
    double d = 2./(h*h); // diagonal constant

    for(int i=0 ; i<n ; i++){
            for(int j=0 ; j<n ; j++){
                if(j==i)
                    A[i][j] = d + V[i]; // diagonal matrix elements
                else if(j == i+1 || j == i-1)
                    A[i][j] = e; // non-diagonal non-zero matrix elements
                else
                    A[i][j] = 0.0;
            }
        }

    cout << A << endl;


    // Doing the calculations
    int iterations = 0;
    double epsilon = 1.0e8; // tolerance
    double max_number_iterations = (double) n * (double) n * (double) n;
    double max_offdiag = maxoffdiag(A, &k, &l, n);

    while (fabs (max_offdiag) > epsilon && (double) iterations < max_number_iterations){
        max_offdiag = maxoffdiag(A, &k, &l, n); // find max. matrix element
        rotate(A, R, k, l, n); // rotate the matrix
        iterations++;
    }

    cout << "Number of iterations: " << iterations << endl;

    finish = clock(); // final time
    cout << "Time: " << "\t" << ((finish - start)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time

    // Clear memory
    for (int i = 0; i < n ; i++)
        delete [] A[i];
    delete [] A;

    return 0;
}
