// Leandro Panlilio Viray
// ECE 3340 
// Numerical Methods for Electrical and Computer Engineers
// Programming Assignment 1
// October 3, 2018

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Initialize Global Variables for use in the functions and main loop
double Up_Bound = 0.0;					// Initialize Upper Bound value for application of methods.
double Low_Bound = 0.0;					// Initialize Lower Bound value for application of methods.
double *func;							// Initialize pointer to polynomial and interval array. 
double t = 0.00001;						// Initialize tolerance value (TOL) 10^-5 , for root value comparison.
double midpoint = 0.0;					// Intialize midpoint value to calculate roots.
double root = 0.0;						// Initialize root estimate of polynomial.
double root_n = 0.0;					// Initialize next root estimate of polynomial.
int N_max = 0;							// Initialize maximum number of iterations for methods.
int order = 0;							// Initialize degree of polynomial count.
double a; 								// Initialize Lower Bound symbol.
double b;								// Initialize Upper Bound symbol. 
double c;								// Initialize Midpoint symbol. 
int bisections;							// Initialize bisection counter. 
int flag;								// Initialize flag.

// Intialize Function Prototypes
double Horner(double nC);				// Initialize Horner's Algorithm Method. 
double Horner_Derivative(double nC);	// Initialize Horner's Derivative Algorithm Method.
double Newton();						// Initialize Newton's Method.
double Bisection();						// Initialize Bisection Method.

// Horner's Algorithm Method 
double Horner(double nC) {

	double x = 0.0;

	for (int i = 0; i < order; i++) {    // We apply synthetic division here to evaluate the value of the polynomial.
		x = x * nC + func[i];
	}

	return x;
}

// Horner's Algorithm Method for Derivatives
double Horner_Derivative(double nC) {

	double y = 0.0;
	double x = 0.0;

	for (int i = 0; i < order; i++) {		//  We apply synthetic divison here to evaluate the value of the derivative of the polynomial. 
		y = y * nC + x;
		x = x * nC + func[i];
	}

	return y;

}

// Newton's Method 
double Newton() {

	if ((Horner(Low_Bound) * Horner(Up_Bound)) > 0) {									// Check if the interval even has a root within it. 
		printf("\nError, there is no root found within this interval.\n\n");
		return 0;
	}

	double R1, R2;
	N_max = (int)ceil((log((double)(Up_Bound - Low_Bound)) - log(t)) / log(2.0));		// Estimate maximum number of iterations to apply for Newton's and Bisection Methods.  
	midpoint = ((Up_Bound + Low_Bound) / 2.0);
	root = midpoint;																	// Initialize root estimate from given interval.

	int i = 0;																			// Initialize iteration counter. 

	if ((root >= Low_Bound) && (root <= Up_Bound)) {									// Make sure roots stay within bounds. 
		for (i = 0; i <= (N_max + 1); i++) {

			if (i > 2) {
				if ((R1 == root) && (R2 == root)) {						// Check if differents roots end up being the same value to know if we need to use the Bisection method instead.
					R1 = root;
					R2 = root;
					root = Bisection();									// Use Bisection method instead.
					i = 0;
					continue;											// Continue on to next iteration after getting new root from Bisection. 
				}
			}


			if ((root < Low_Bound || root > Up_Bound) || fabs(Horner_Derivative(root)) == 0) { 		// Check if the derivative of the root is equal to 0 or the roots are out of bounds to know if we need to use the Bisection method instead.			     
				R1 = root;
				R2 = root;
				root = Bisection();										// Use Bisection method instead. 
				i = 0;
				continue;												// Continue on to next iteration after getting new root from Bisection. 
			}

			else {
				R2 = R1;
				R1 = root;
				root = root - (Horner(root) / Horner_Derivative(root));					// Implement Newton's Method to solve for root of polynomial.

			}


			if (fabs(Horner(root)) < t) {												// Break out if found root early before reaching maximum iterations.
				break;
			}

			if (i == (N_max + 1)) {														// If Newton's Method exceeds maximum iterations, use Bisection Method instead. 
				R1 = root;
				R2 = root;
				root = Bisection();														// Use Bisection method instead. 
				i = 0;
				continue;																// Continue on to next iteration after getting new root from Bisection.
			}
		}
	}


	printf("\nThe root of the polynomial utilizing the Newton's Method is %f with an iteration number of %d. %d bisections were applied.\n\n", root, i, bisections);		// Announce final root estimation regardless of hitting the maximum number of iterations or not.

	return root;																		// Return final root back into main function.

}

// Bisection Method
double Bisection() {


	c = ((a + b) / 2.0);													// Use midpoint formula on interval [a,b] to determine root estimate.
	bisections++;

	if ((Horner(c) * Horner(a)) > 0.0) {									// Determine if new root estimate should replace upper or lower bound depending on f(c) sign value. 

		a = c;

	}

	else {

		b = c;																// New interval is created here for use in next iteration of Bisection Method. 

	}


	return c;

}



int main(int argc, char** argv) {

	int n = 1;												// Set up starting point of array, past address. 
	order = argc - 3;										// Set up maximum elements for storage of coefficients of polynomial.
	int flag;												// Initialize flag to know what part of argument the program needs to store and print to screen on Unix server. 

	func = (double*)malloc(order * sizeof(double));			// Allocate storage to take in the polynomial function.

	while (n < argc) {

		if (n < order) {									// Flag 1 stores the coefficients of polynomial function.

			flag = 1;

		}

		if (n == order) {

			flag = 2;										// Flag 2 stores the constant of polynomial function. 

		}

		if (n == (order + 1)) {

			flag = 3;										// Flag 3 stores the lower bound of interval of polynomial function. 

		}

		if (n == (order + 2)) {

			flag = 4;										// Flag 4 stores the upper bound of interval of polynomial function. 

		}

		switch (flag) {
		case 1:
			func[n - 1] = atof(argv[n]);
			printf("%fx^%d + ", func[n - 1], (order - n));	// Print inputted polynomial function to screen.
			break;
		case 2:
			func[n - 1] = atof(argv[n]);
			printf("%f\n", func[n - 1]);					// Print constant of polynomial function to screen. 
			break;
		case 3:
			Low_Bound = atof(argv[n]);						// Store lower bound of inputted interval. 
			break;
		case 4:
			Up_Bound = atof(argv[n]);						// Store upper bound of inputted interval. 
			break;
		}
		n++;
	}

	a = Low_Bound;											// Making sure lower bound initialized properly in main loop.
	b = Up_Bound;											// Making sure upper bound initialized properly in main loop. 


	N_max = (int)ceil((log((double)(Up_Bound - Low_Bound)) - log(t)) / log(2.0));		// Estimate maximum number of iterations to apply for Newton's and Bisection Methods.  

	if (midpoint > Up_Bound && midpoint < Low_Bound) {
		printf("\n%lf is an invalid root because it is out of the determined interval.\n", midpoint);	// Set up case in which the estimated root is out of bounds, invalid estimate.  
	}

	else {

		Newton();																		// Apply Newtons or Bisection Method to estimate root of polynomial function. 

	}


	return 0;

}


