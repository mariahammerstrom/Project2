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
#include <algorithm>
#include <armadillo>
#include "lib.h"

using namespace std;
using namespace arma;



// FUNCTIONS

// Calculates and return the value of the potential for a given value of x.
double potential_HO(double x){
    return x*x;
}


double potential_C(double x, double omega){
    return omega*omega*x*x + 1./x;
}


// Find maximum matrix element
double maxoffdiag(mat A, int * k, int * l, int n){
    double max = 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            if (fabs(A(i,j)) > max){
                max = fabs(A(i,j)); // value of max element
                *l = i; // row index of max element
                *k = j; // column index of max element
            }
        }
    }
    return max;
}

// Find values of cos and sin
mat rotate (mat A, mat R, int k, int l, int n){
    double s,c; // s = sin(theta), c = cos(theta)
    if (A(k,l) != 0.0){
        double t, tau; // t = tan(theta), tau = cot(2*theta)
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/(- tau + sqrt(1.0 + tau*tau));
        }
        c = 1./sqrt(1 + t*t);
        s = c*t;
    } else {
        c = 1.0;
        s = 0.0;
    }

    // Changing the matrix elements with indices k and l
    double a_kk = A(k,k);
    double a_ll = A(l,l);

    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    // Changing the remaining elements
    for (int i = 0; i < n; i++){
        if (i != k && i != l){
            double a_ik = A(i,k);
            double a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }

        // Compute new eigenvectors
        double r_ik = R(i,k);
        double r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return A;
}




// MAIN PROGRAM
int main()
{
    // CONSTANTS
    int n = 200; // number of steps
    double rho_min = 0.0; // minimum value
    double rho_max = 10.0; // maximum value (set by the user)
    double h = (rho_max - rho_min)/n; // step length
    double omega = 0.01;

    double e = -1.0/(h*h); // non-diagonal constant
    double d = 2.0/(h*h); // diagonal constant

    cout << "INITIAL CONDITIONS:" << endl;
    cout << "rho_min = " << rho_min << endl;
    cout << "rho_max = " << rho_max << endl;
    cout << "No. of steps = " << n << endl;
    cout << "Omega = " << omega << endl;
    cout << endl;


    // ARMADILLO SOLVER
    cout << "ARMADILLO:" << endl;
    clock_t start_armadillo, finish_armadillo; // declare start and final time for Armadillo solver
    start_armadillo = clock();

    // Set up dimensionless variable rho and potential
    vec rho(n+1);
    vec V(n+1);

    for (int i = 0; i < n+1; i++){
        rho(i) = rho_min + (i)*h;
        V(i) = potential_C(rho(i),omega);
    }

    // Set up tridiagonal matrix A
    mat A = zeros<mat>(n,n);

    A(0,0) = d + V(1); // upper-left element
    A(0,1) = e;

    for (int i = 1; i < n-1; i++){
        A(i,i-1) = e; // non-diagonal non-zero matrix element
        A(i,i) = d + V(i+1); // diagonal elements
        A(i,i+1) = e; // non-diagonal non-zero matrix element
    }

    A(n-1,n-2) = e;
    A(n-1,n-1) = d + V(n);
    mat B = A; // copy matrix for use with our algorithm


    // EIGENVALUES
    cout << "Lowest eigenvalues:" << endl;
    vec eigval;
    mat eigvec;

    eig_sym(eigval,eigvec,A);

    double omega_e = sqrt(3)*omega;
    double V_0 = 2*((3./2)*pow(omega/2.,2./3));

    for (int i = 0; i < 3; i++){
        double eigval_analytic = V_0 + omega_e*(i + 0.5); // analytical eigenvalue
        cout << "Numeric: " << eigval(i) << "\t" << "Analytic: " << eigval_analytic << endl;
    }

    finish_armadillo = clock(); // final time
    cout << "Time: " << "\t" << ((finish_armadillo - start_armadillo)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time


    // 1st EIGENVECTOR
    ofstream file("Eigenvectors_" + to_string(omega) + ".txt"); // File for wave function
    ofstream file2("Eigenvectors2_" + to_string(omega) + ".txt"); // File for wave function squared

    for(int i = 0; i<n; i++){
        file << rho(i) << "\t" << eigvec(i,0) << endl; // Write solution to file
        file2 << rho(i) << "\t" << eigvec(i,0)*eigvec(i,0) << endl; // Write solution to file
    }

    file.close();
    file2.close();


    // OUR ALGORITHM
    cout << endl << "OUR ALGORITHM:" << endl;
    clock_t start, finish; // declare start and final time
    start = clock();

    // Set up eigenvector matrix
    mat R = zeros<mat>(n,n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                R(i,j) = 1.0; // diagonal matrix elements
            } else {
                R(i,j) = 0.0; // non-diagonal matrix elements
            }
        }
    }

    int iterations = 0;
    int k,l;
    double epsilon = 1.0e-10; // tolerance
    double max_number_iterations = 50 * (double) n * (double) n * (double) n;
    cout << "Max no. iterations: " << max_number_iterations << endl;

    double max_offdiag = maxoffdiag(B, &k, &l, n);

    while (fabs (max_offdiag) > epsilon && (double) iterations < max_number_iterations){
        B = rotate(B, R, k, l, n); // rotate the matrix
        max_offdiag = maxoffdiag(B, &k, &l, n); // find max. matrix element
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;


    // Eigenvalues
    vec Eigval(n);
    for (int i = 0; i < n; i++){
        Eigval(i) = B(i,i);
    }

    cout << "Eigenvalues:" << endl;
    Eigval = sort(Eigval,"ascend");
    for (int i = 0; i < 3; i++){
        cout << Eigval(i) << endl; // Print out first 3 eigenvalues
    }

    finish = clock(); // final time
    cout << "Time: " << "\t" << ((finish - start)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time

    return 0;
}
