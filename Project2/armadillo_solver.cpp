// Using Armadillo to solve an eigenvalue problem
// Source: https://github.com/CompPhysics/ComputationalPhysics1/blob/gh-pages/doc/pub/eigvalues/pdf/eigvalues-print.pdf

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
double potential_HO(double x){
    return x*x; // harmonic oscillator potential
}

// Calculates and return the value of the potential for a given value of x and omega.
double potential_coulomb(double x, double omega){
    return omega*omega*x*x + 1./x; // harmonic oscillator potential with coulomb interaction
}



// MAIN PROGRAM
int main()
{
    clock_t start, finish; // declare start and final time
    start = clock();

    int n = 10; // number of steps
    double rho_min = 0.0; // minimum value
    double rho_max = 10.0; // maximum value (set by the user)
    double h = (rho_max - rho_min)/( (double) n); // step length


    // Set up potential
    vec rho(n);
    vec V(n);

    for (int i = 0; i < n; i++){
        rho(i) = rho_min + i*h;
        V(i) = potential_HO(rho(i));
    }


    // Set up tridiagonal matrix A
    mat A = zeros<mat>(n,n);

    double e = -1./(h*h); // non-diagonal constant
    double d = 2./(h*h); // diagonal constant

    A(0,0) = 2./(h*h) + V(0); // upper-left element
    A(0,1) = e;

    for (int i = 1; i < n; i++){
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
    cout << "rho_min = " << rho_min << endl;
    cout << "rho_max = " << rho_max << endl;
    cout << "No. of steps = " << n << endl;
    cout << "Five lowest eigenvalues:" << endl;
    for (int i = 0; i < 5; i++){
        cout << eigval(i) << endl; // correct thing to print out??
    }

    finish = clock(); // final time
    cout << "Time: " << "\t" << ((finish - start)/CLOCKS_PER_SEC) << " seconds" << endl; // print elapsed time


    return 0;
}
