/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 #include <stdio.h>
 #include <ctype.h>
 #include <stdlib.h>
 #include <time.h>
 #include <math.h>
 #include "functions_data.h"
 #include "functions_prng.h"

#include "functions.h"

#ifndef M_PI
	#define M_PI acos(-1.0)
#endif

#ifndef M_E
	#define M_E 2.7182818284590452353602874713527
#endif


static const double griewank_m_d2[2][2] = GRIEWANK_M_D2;
static const double griewank_m_d10[10][10] = GRIEWANK_M_D10;
static const double griewank_m_d30[30][30] = GRIEWANK_M_D30;
static const double griewank_m_d50[50][50] = GRIEWANK_M_D50;

/// Sphere function
double Sphere(double * vars, unsigned int nvars) {

	unsigned int i;
	double fitness = 0.0;

	for (i = 0; i < nvars; ++i) {
		fitness += (double) (vars[i] * vars[i]);
	}

	return fitness;
}

/// Quadric function
double Quadric(double * vars, unsigned int nvars) {

	unsigned int i, j;
	double fitaux, fitness = 0.0;

	for (i = 1; i <= nvars; ++i) {

		fitaux = 0;

		for (j = 0; j < i; ++j) {
			fitaux += (double) vars[j];
		}

		fitness += fitaux * fitaux;
	}
	return fitness;
}

/// Hyperellipsoid function
double Hyper(double * vars, unsigned int nvars) {

	unsigned int i;
	double fitness = 0.0;

	for (i = 0; i < nvars; ++i) {
		fitness += i * (double) (vars[i] * vars[i]);
	}

	return fitness;
}

/// Rastringin function
double Rastrigin(double * vars, unsigned int nvars) {

	unsigned int i;
	double fitness = 0.0;

	for (i = 0; i < nvars; ++i) {
		fitness += (vars[i] * vars[i] - 10 *
			cos(2.0 * M_PI * vars[i]) + 10);
	}

	return fitness;
}

/// Griewank function
double Griewank(double * vars, unsigned int nvars) {

	unsigned int i;
	double fitness1, fitness2, fitness;

	fitness1 = 0.0;
	fitness2 = 1.0;

	for (i = 0; i < nvars; ++i) {
		fitness1 += (double) (vars[i] * vars[i]);
		fitness2 *= cos(vars[i] / sqrt(i + 1.0));
	}

	fitness = 1 + (fitness1 / 4000) - fitness2;
	return fitness;
}

/// Schaffer function No. 6
double Schaffer6(double * vars, unsigned int nvars) {

	double x, y;
	double temp1, temp2;
	(void) nvars;

	x = (double) vars[0];
	y = (double) vars[1];

	temp1 = (double) sin(sqrt(x * x + y * y));
	temp2 = (double) (1 + 0.001 * (x * x + y * y));

	return (double) (0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2));
}

/// Weierstrass function
double Weierstrass(double * vars, unsigned int nvars) {

	unsigned int i, j;
	double res;
	double sum;
	double a, b;
	unsigned int k_max;

	a = 0.5;
	b = 3.0;
	k_max = 20;
	res = 0.0;

	for (i = 0; i < nvars; ++i) {
		sum = 0.0;
		for (j = 0; j <= k_max; j++)
			sum += pow(a, j) *
				cos(2.0 * M_PI * pow(b, j) *
				(vars[i] + 0.5));
		res += sum;
	}

	sum = 0.0;
	for (j = 0; j <= k_max; ++j) {
		sum += pow(a,j) * cos(2.0 * M_PI * pow(b, j) * 0.5);
	}

	return res - nvars * sum;
}

/// Ackley function
double Ackley(double * vars, unsigned int nvars) {

	double fitness;
	unsigned int j;
	double fitaux1, fitaux2;

	fitness = 0.0;
	fitaux1 = 0;
	fitaux2 = 0;

	for (j = 0; j < nvars; ++j) {
		fitaux1 += vars[j] * vars[j];
		fitaux2 += cos(2 * M_PI * vars[j]);
	}

	fitness = -20 * exp(-0.2 * sqrt(fitaux1 / nvars))
		- exp(fitaux2 / nvars) + 20 + M_E;

	return fitness;
}

/// Shifted Quadric With Noise function, also know as Shifted Schwefel 2 w/ Noise
double ShiftedQuadricWithNoise(double * vars, unsigned int nvars) {

	const double o[] = SCHWEFEL_102_DATA;
	unsigned int i, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	float shifted_vars[nvars];

	for (i = 0; i < nvars; ++i) {
		shifted_vars[i] = vars[i] - o[i];
	}

	for (i = 0; i < nvars; ++i) {
		sum2 = 0.0;
		for (j = 0; j <= i; ++j) {
			sum2 += shifted_vars[j];
		}
		sum1 += sum2 * sum2;
	}
	sum1 = sum1 * (1.0 + 0.4 * fabs(runif01(vars, nvars)));
	return sum1;
}

/// Rotated Griewank function
double RotatedGriewank(double * vars, unsigned int nvars) {

	const double (*m)[nvars];
	unsigned int i, j;
	double rotated_vars[nvars];

	switch (nvars) {
		case 2: m = griewank_m_d2; break;
		case 10: m = griewank_m_d10; break;
		case 30: m = griewank_m_d30; break;
		case 50: m = griewank_m_d50; break;
		default:
			fprintf(stderr, "The specified number of variables/dimensions (%u)"
				"is not supported by the '%s' function.\n", nvars, __func__);
			return 0.0;
	}

	for (i = 0; i < nvars; ++i) {
		rotated_vars[i] = 0.0;
		for (j = 0; j < nvars; j++) {
			rotated_vars[i] += (*(m + i))[j] * vars[i];
		}
	}

	return Griewank(rotated_vars, nvars);
}

// Basic function used to test the internal noise-generating PRNG
double RandomTest(double * vars, unsigned int nvars) {
	return runif01(vars, nvars);
}
