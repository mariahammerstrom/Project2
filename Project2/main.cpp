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
    //omega = 0.01;
    //return omega*omega*x*x + 1./x;
    return x*x;
}

// Find maximum matrix element
double maxoffdiag(double ** A, int * k, int * l, int n){
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
void rotate (double ** A, double ** R, int k, int l, int n){
    double s,c; // s = sin(theta), c = cos(theta)
    if (A[k][l] != 0.0){
        double t, tau; // t = tan(theta), tau = cot(2*theta)
        tau = (A[l][l] - A[k][k]/(2*A[k][l]));
        if (tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = - 1.0/(- tau + sqrt(1.0 + tau*tau));
        }
        c = 1./sqrt(1 + t*t);
        s = c*t;
    } else {
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

void jacobi_method(double ** A, double ** R, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                R[i][j] = 1.0; // diagonal matrix elements
            } else {
                R[i][j] = 0.0; // non-diagonal matrix elements
            }
        }
    }

    int iterations = 0;
    int k,l;
    double epsilon = 1.0e-8; // tolerance
    double max_number_iterations = (double) n * (double) n * (double) n;
    double max_offdiag = maxoffdiag(A, &k, &l, n);

    while (fabs (max_offdiag) > epsilon && (double) iterations < max_number_iterations){
        max_offdiag = maxoffdiag(A, &k, &l, n); // find max. matrix element
        rotate(A, R, k, l, n); // rotate the matrix
        cout << "Rotation #" << iterations << ":" << endl;
        cout << "B[0][0] = " << A[0][0] << endl;
        cout << "B[0][1] = " << A[0][1] << endl;
        cout << "B[1][0] = " << A[1][0] << endl;
        cout << "B[1][1] = " << A[1][1] << endl;
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;
    return;
}



// MAIN PROGRAM
int main()
{
    // CONSTANTS
    int n = 2; // number of steps
    double rho_min = 0.0; // minimum value
    double rho_max = 3.0; // maximum value (set by the user)
    double h = (rho_max - rho_min)/( (double) n); // step length

    double e = -1./(h*h); // non-diagonal constant
    double d = 2./(h*h); // diagonal constant

    cout << "rho_min = " << rho_min << endl;
    cout << "rho_max = " << rho_max << endl;
    cout << "No. of steps = " << n << endl;
    cout << endl;


    // ARMADILLO SOLVER
    cout << "ARMADILLO" << endl;
    clock_t start_armadillo, finish_armadillo; // declare start and final time
    start_armadillo = clock();

    // Set up potential
    vec rho(n);
    vec V(n);

    for (int i = 0; i < n; i++){
        rho(i) = rho_min + i*h;
        V(i) = potential(rho(i));
    }

    // Set up tridiagonal matrix A
    mat A = zeros<mat>(n,n);

    A(0,0) = 2./(h*h) + V(0); // upper-left element
    A(0,1) = e;

    for (int i = 1; i < n-1; i++){
        rho = rho_min + i*h;
        A(i,i-1) = e; // non-diagonal non-zero matrix element
        A(i,i) = d + V(i); // diagonal elements
        A(i,i+1) = e; // non-diagonal non-zero matrix element
    }

    A(n-1,n-2) = e;
    A(n-1,n-1) = d + V(n-1);

    // Diagonalize and find eigenvalues
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);

    // Print out results
    cout << "Lowest eigenvalues:" << endl;
    for (int i = 0; i < 2; i++){
        cout << eigval(i) << endl; // correct thing to print out??
    }

    finish_armadillo = clock(); // final time
    cout << "Time: " << "\t" << ((finish_armadillo - start_armadillo)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time



    // OUR ALGORITHM
    cout << endl << "OUR ALGORITHM" << endl;
    clock_t start, finish; // declare start and final time
    start = clock();

    // Set up potential
    double rho_i,V_i[n];

    for (int i = 0; i < n; i++){
        rho_i = rho_min + i*h;
        V_i[i] = potential(rho_i);
    }

    // Set up eigenvector matrix
    double ** R;
    R = new double*[n];
    for (int i = 0; i < n; i++)
        R[i] = new double[n];

    // Set up tridiagonal matrix
    double ** B;
    B = new double*[n];
    for (int i = 0; i < n; i++)
        B[i] = new double[n];

    for(int i=0 ; i<n ; i++){
            for(int j=0 ; j<n ; j++){
                if(j==i)
                    B[i][j] = d + V_i[i]; // diagonal matrix elements
                else if(j == i+1 || j == i-1)
                    B[i][j] = e; // non-diagonal non-zero matrix elements
                else
                    B[i][j] = 0.0;
            }
        }
    cout << "Before Jacobi:" << endl;
    cout << "B[0][0] = " << B[0][0] << endl;
    cout << "B[0][1] = " << B[0][1] << endl;
    cout << "B[1][0] = " << B[1][0] << endl;
    cout << "B[1][1] = " << B[1][1] << endl;

    jacobi_method(B,R,n);
    cout << "Lowest eigenvalues:" << endl;
    for (int i=0; i<2 ; i++){
        cout << B[i][i] << endl;
    }

    //cout << "After Jacobi:" << endl;
    //cout << "B[0][0] = " << B[0][0] << endl;
    //cout << "B[0][1] = " << B[0][1] << endl;
    //cout << "B[1][0] = " << B[1][0] << endl;
    //cout << "B[1][1] = " << B[1][1] << endl;

    // Clear memory
    for (int i = 0; i < n ; i++)
        delete [] B[i];
    delete [] B;

    finish = clock(); // final time
    cout << "Time: " << "\t" << ((finish - start)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time

    return 0;
}
