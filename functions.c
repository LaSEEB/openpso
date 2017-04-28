/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "functions.h"
#include <math.h>

#ifndef M_PI
	#define M_PI acos(-1.0)
#endif

#ifndef M_E
	#define M_E 2.7182818284590452353602874713527
#endif

/// Sphere function
long double Sphere(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i;
	long double fitness = 0.0;
	(void) prng_state;

	for (i = 0; i < nvars; ++i) {
		fitness += (long double) (position[i] * position[i]);
	}

	return fitness;
}

/// Quadric function
long double Quadric(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i, j;
	long double fitaux, fitness = 0.0;
	(void) prng_state;

	for (i = 1; i <= nvars; ++i) {

		fitaux = 0;

		for (j = 0; j < i; ++j) {
			fitaux += (long double) position[j];
		}

		fitness += fitaux * fitaux;
	}
	return fitness;
}

/// Hyperellipsoid function
long double Hyper(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i;
	long double fitness = 0.0;
	(void) prng_state;

	for (i = 0; i < nvars; ++i) {
		fitness += i * (long double) (position[i] * position[i]);
	}

	return fitness;
}

/// Rastringin function
long double Rastringin(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i;
	long double fitness = 0.0;
	(void) prng_state;

	for (i = 0; i < nvars; ++i) {
		fitness += (position[i] * position[i] - 10 *
			cos(2.0 * M_PI * position[i]) + 10);
	}

	return fitness;
}

/// Griewank function
long double Griewank(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i;
	long double fitness1, fitness2, fitness;
	(void) prng_state;

	fitness1 = 0.0;
	fitness2 = 1.0;

	for (i = 1; i < nvars; ++i) {
		fitness1 += (long double) (position[i] * position[i]);
		fitness2 *= cos(position[i] / sqrt(i + 1.0));
	}

	fitness = 1 + (fitness1 / 4000) - fitness2;
	return fitness;
}

/// Schaffer function
long double Schaffer(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	long double x, y;
	long double temp1, temp2;
	(void) prng_state;
	(void) nvars;

	x = (long double) position[0];
	y = (long double) position[1];

	temp1 = (long double) sin(sqrt(x * x + y * y));
	temp2 = (long double) (1 + 0.001 * (x * x + y * y));

	return (long double) (0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2));
}

/// Weierstrass function
long double Weierstrass(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	unsigned int i, j;
	long double res;
	long double sum;
	long double a, b;
	unsigned int k_max;
	(void) prng_state;

	a = 0.5;
	b = 3.0;
	k_max = 20;
	res = 0.0;

	for (i = 0; i < nvars; ++i) {
		sum = 0.0;
		for (j = 0; j <= k_max; j++)
			sum += pow(a, j) *
				cos(2.0 * M_PI * pow(b, j) *
				(position[i] + 0.5));
		res += sum;
	}

	sum = 0.0;
	for (j = 0; j <= k_max; ++j) {
		sum += pow(a,j) * cos(2.0 * M_PI * pow(b, j) * 0.5);
	}

	return res - nvars * sum;
}

/// Ackley function
long double Ackley(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	long double fitness;
	unsigned int j;
	long double fitaux1, fitaux2;
	(void) prng_state;

	fitness = 0.0;
	fitaux1 = 0;
	fitaux2 = 0;

	for (j = 0; j < nvars; ++j) {
		fitaux1 += position[j] * position[j];
		fitaux2 += cos(2 * M_PI * position[j]);
	}

	fitness = -20 * exp(-0.2 * sqrt(fitaux1 / nvars))
		- exp(fitaux2 / nvars) + 20 + M_E;

	return fitness;
}

/// Shifted Quadric With Noise function, also know as Shifted Schwefel 2 w/ Noise
long double ShiftedQuadricWithNoise(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	const long double o[] = SCHWEFEL_102_DATA;
	unsigned int i, j;
	long double sum1 = 0.0;
	long double sum2 = 0.0;
	float shifted_position[nvars];

	for (i = 0; i < nvars; ++i) {
		shifted_position[i] = position[i] - o[i];
	}

	for (i = 0; i < nvars; ++i) {
		sum2 = 0.0;
		for (j = 0; j <= i; ++j) {
			sum2 += shifted_position[j];
		}
		sum1 += sum2 * sum2;
	}
	sum1 = sum1 * (1.0 + 0.4 * fabs(rds_luniform(prng_state, 0.0, 1.0)));
	return sum1;
}

/// Rotated Griewank function
long double RotatedGriewank(PARTICLE position, unsigned int nvars, mt_state * prng_state) {

	const long double m[50][50] = GRIEWANK_M_D50;
	unsigned int i, j;
	long double fitness1, fitness2, fitness;
	float rotated_position[nvars];
	(void) prng_state;

	for (j = 0; j < nvars; ++j) {
		rotated_position[j] = 0.0;
		for (i = 0; i < nvars; i++) {
			rotated_position[j] += m[i][j] * position[i];
		}
	}

	fitness1 = 0.0;
	fitness2 = 1.0;
	for (i = 1; i < nvars + 1; ++i) {
		fitness1 = fitness1 + (long double) (position[i] * position[i]);
		fitness2 = fitness2 * (cos(position[i] / sqrt((float)i)));
	}
	fitness = 1 + (fitness1 / 4000) - fitness2;
	return fitness;
}


SelFunc getSelFunc(unsigned int problem) {
	switch (problem) {
		case 1: return &Sphere;
		case 2: return &Quadric;
		case 3: return &Hyper;
		case 4: return &Rastringin;
		case 5: return &Griewank;
		case 6: return &Schaffer;
		case 7: return &Weierstrass;
		case 8: return &Ackley;
		case 9: return &ShiftedQuadricWithNoise;
		case 10: return &RotatedGriewank;
		default: ERROR_EXIT("Unknown function!\n");
	}
}
