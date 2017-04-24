/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "functions.h"
#include <math.h>

#ifndef M_PI
	#define M_PI acos(-1.0)
#endif

// SPHERE FUNCTION
long double Sphere(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness = 0.0;
	for (i = 0; i < nvars; ++i)
		fitness = fitness + (long double) (position[i] * position[i]);
	return fitness;
}

// HYPERELIPSOID FUNCTION
long double Hyper(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness = 0.0;
	for (i = 0; i < nvars; ++i)
		fitness = fitness + i * (long double) (position[i] * position[i]);
	return fitness;
}

// HYPERELIPSOID FUNCTION
long double Quartic(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness = 0.0;
	for (i = 0; i < nvars; ++i)
		fitness = fitness + i * (long double)
		(position[i] * position[i] * position[i] * position[i]);
	return fitness;
}

// ACKLEY FUNCTION
long double Ackley(PARTICLE position, unsigned int nvars) {
	long double fitness;
	unsigned int j;
	long double fitaux1, fitaux2;
	fitness = 0.0;
	fitaux1 = 0;
	fitaux2 = 0;
	for (j = 0; j < nvars; ++j) {
		fitaux1 = fitaux1 + position[j] * position[j];
		fitaux2 = fitaux2 + cos(2 * M_PI * position[j]);
	}
	fitness = (22.71828182846 - 20 * pow(M_PI,
		(-0.2 * (sqrt((fitaux1 / (long double) nvars))))) - pow(2.71828182846,
		(fitaux2 / (long double) nvars)));
	return fitness;
}

// ROSENBROCK FUNCTION
long double Rosenbrock(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness, tt, x, xx;
	fitness = 0.0;
	for (i = 0; i < nvars - 1; ++i) {
		tt = position[i] * position[i];
		x = position[i];
		xx = position[i + 1];
		fitness = fitness + 100 * (xx - tt) * (xx - tt) + (1 - x) * (1 - x);
	}
	return fitness;
}

// RASTRIGIN FUNCTION
long double Rastringin(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness = 0.0;
	for (i = 0; i < nvars; ++i)
		fitness = fitness +
			(position[i] * position[i] - 10 *
			cos(2.0 * M_PI * position[i]) + 10);
	return fitness;
}

long double Griewank(PARTICLE position, unsigned int nvars) {
	unsigned int i;
	long double fitness1, fitness2, fitness;
	fitness1 = 0.0;
	fitness2 = 1.0;
	for (i = 1; i < nvars + 1; ++i) {
		fitness1 = fitness1 + (long double) (position[i] * position[i]);
		fitness2 = fitness2 * (cos(position[i] / sqrt((float) i)));
	}
	fitness = 1 + (fitness1 / 4000) - fitness2;
	return fitness;
}

long double Schaffer(PARTICLE position, unsigned int nvars) {
	long double x, y;
	long double temp1, temp2;
	(void) nvars;
	x = (long double) position[0];
	y = (long double) position[1];
	temp1 = (long double) sin(sqrt(x * x + y * y));
	temp2 = (long double) (1 + 0.001 * (x * x + y * y));
	return (long double) (0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2));
}

long double Weierstrass (PARTICLE position, unsigned int nvars) {
	unsigned int i, j;
	long double res;
	long double sum;
	long double a, b;
	unsigned int k_max;
	a = 0.5;
	b = 3.0;
	k_max = 20;
	res = 0.0;
	for (i = 0; i < nvars; i++) {
		sum = 0.0;
		for (j = 0; j <= k_max; j++)
			sum += pow(a, j) * cos(2.0 * M_PI * pow(b,j) * (position[i] + 0.5));
		res += sum;
	}
	sum = 0;
	//for (j = 0; j <= k_max; j++)
	//	sum = sum+pow(a,j)*cos(2.0*M_PI*pow(b,j)*0.5);
	res = res + 60;
	return (res);
}

long double Schwefel (PARTICLE position, unsigned int nvars) {
	unsigned int j;
	long double fitness;
	long double sum;
	sum = 0.0;
	for (j = 0; j < nvars; j++)
		sum = sum + position[j] * sin(sqrt(fabs((float) position[j])));
	fitness = 418.9829 * nvars - sum;
	return (fitness);
}

SelFunc getSelFunc(unsigned int problem) {
	switch (problem) {
		case 1: return &Sphere;
		case 2: return &Rosenbrock;
		case 3: return &Rastringin;
		case 4: return &Griewank;
		case 5: return &Schaffer;
		case 6: return &Weierstrass;
		case 7: return &Ackley;
		case 8: return &Schwefel;
		case 9: return &Hyper;
		case 10: return &Quartic;
		default: ERROR_EXIT("Unknown function!\n");
	}
}
