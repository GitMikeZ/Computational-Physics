// Name: Minqi Zhang
// Student #: 15383078
// Code was revised from Computational Physics's
// exponential derivative by Morten Hjorth-Jensen
// under the Creative Commons License.
//
// Question #1

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void initialise(double*, double*, int*);
void first_deriviative_f(int , double , double , double *, double*);
void first_deriviative_c(int , double , double , double *, double*);
void first_deriviative_e(int , double , double , double *, double*);

void output(char *, double *, double *, double, int);

int main()
{
	// declarations of variables
	int number_of_steps;
	double x, ini_step;
	double *h_step;
	double *computed_derivative_f, *computed_derivative_c, *computed_derivative_e;

	// read input data from screen
	initialise (&ini_step, &x, &number_of_steps);

	// allocate memory for one-dimensional array
	// h_step and computed_derivative
	h_step = malloc(number_of_steps * sizeof(double));

	computed_derivative_f = malloc(number_of_steps*sizeof(double));
	computed_derivative_c = malloc(number_of_steps*sizeof(double));
	computed_derivative_e = malloc(number_of_steps*sizeof(double));

	// compute first derivative of 1/cosh(x) with forward difference
	first_deriviative_f(number_of_steps, x, ini_step, h_step, computed_derivative_f);

	// print results
	output("out_forward.dat", h_step, computed_derivative_f, x, number_of_steps );

	// compute first derivative of 1/cosh(x) with central difference
	first_deriviative_c(number_of_steps, x, ini_step, h_step, computed_derivative_c);

	// print results
	output("out_central.dat", h_step, computed_derivative_c, x, number_of_steps );

	// compute first derivative of 1/cosh(x) with extrapolated difference
	first_deriviative_e(number_of_steps, x, ini_step, h_step, computed_derivative_e);

	// print results
	output("out_extrapolated.dat", h_step, computed_derivative_e, x, number_of_steps );

	// free memory
	free(h_step);
	free(computed_derivative_f);
	free(computed_derivative_c);
	free(computed_derivative_e);

	return 0;
}

void initialise (double *initial_step, double *x, int *number_of_steps)
{
	printf("Read in the initial step, x, and number of steps\n");
	scanf("%lf %lf %d", initial_step, x, number_of_steps);
	return;
} // end of function initalise

void first_deriviative_f( int number_of_steps, double x,
	double ini_step, double *h_step, double *computed_derivive)
{
	int count;
	double h;

	h = ini_step;

	for (count = 0; count < number_of_steps; count++ )
	{
		h_step[count] = h;
		computed_derivive[count] = ((1/cosh(x+h)) - 1/cosh(x))/h;
		h = h*0.5;
	}
} // end of function first_derivative_f

void first_deriviative_c( int number_of_steps, double x,
	double ini_step, double *h_step, double *computed_derivive)
{
	int count;
	double h;

	h = ini_step;

	for (count = 0; count < number_of_steps; count++ )
	{
		h_step[count] = h;
		computed_derivive[count] = ((1/cosh(x+h)) - 1/cosh(x-h))/(2*h);
		h = h*0.5;
	}
} // end of function first_derivative_c

void first_deriviative_e( int number_of_steps, double x,
	double ini_step, double *h_step, double *computed_derivive)
{
	int count;
	double h;

	h = ini_step;

	for (count = 0; count < number_of_steps; count++ )
	{
		h_step[count] = h;
		computed_derivive[count] = ((1/cosh(x-(2*h)) - (8/cosh(x-h)) +
					   (8/cosh(x+h)) - 1/cosh(x+ (2*h))))/(12*h);
		h = h*0.5;
	}
} // end of function first_derivative_e

void output(char *outfilename, double *h_step, double *computed_derivative, double x, int number_of_steps)
{
	int i;
	FILE *output_file;
	output_file = fopen(outfilename, "w");
	for ( i=0; i < number_of_steps; i++ )
	{
		fprintf(output_file, "%12.5E %12.10E %12.10E %12.5E %12.5E \n",
		h_step[i], computed_derivative[i], (1/cosh(x)),
		log10(h_step[i]),
		log10(fabs((computed_derivative[i]+(sinh(x)/(cosh(x)*cosh(x))))/(-1*sinh(x)/(cosh(x)*cosh(x))))));
	}

	fclose(output_file);
} // end of function output