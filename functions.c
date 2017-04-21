/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "functions.h"
#include <math.h>

extern int numberVariables;
extern int problem;

// SPHERE FUNCTION
long double Sphere (PARTICLE position) {
	int i;
	long double fitness = 0.0;
	for (i = 0;    i < numberVariables;      ++i)
	   fitness = fitness+(long double)(position[i]*position[i]);
	return fitness;
}

// HYPERELIPSOID FUNCTION
long double Hyper (PARTICLE position) {
	int i;
	long double fitness = 0.0;
	for (i = 0;    i < numberVariables;      ++i)
	   fitness = fitness+i*(long double)(position[i]*position[i]);
	return fitness;
}

// HYPERELIPSOID FUNCTION
long double Quartic (PARTICLE position) {
	int i;
	long double fitness = 0.0;
	for (i = 0;    i < numberVariables;      ++i)
	   fitness = fitness+i*(long double)(position[i]*position[i]*position[i]*position[i]);
	return fitness;
}

// ACKLEY FUNCTION
long double Ackley (PARTICLE position) {
	long double fitness;
	int j;
	long double fitaux1, fitaux2;
	long double pi = 3.1415926535;
	fitness = 0.0;
	fitaux1 = 0;
	fitaux2 = 0;
	for (j = 0;	j < numberVariables;	++j) {
		fitaux1 = fitaux1+position[j]*position[j];
		fitaux2 = fitaux2+cos(2*pi*position[j]);
	}
	fitness = (22.71828182846-20*pow(pi, (-0.2*(sqrt((fitaux1/(long double)numberVariables)))))-pow(2.71828182846,(fitaux2/(long double)numberVariables)));
	return fitness;
}

// ROSENBROCK FUNCTION
long double Rosenbrock (PARTICLE position) {
	int i;
	long double fitness, tt, x, xx;
	fitness = 0.0;
	for (i = 0;    i < numberVariables-1;    ++i) {
		tt = position[i]*position[i];
		x = position[i];
		xx = position[i+1];
		fitness = fitness+100*(xx-tt)*(xx-tt)+(1-x)*(1-x);
	}
	return fitness;
}

// RASTRIGIN FUNCTION
long double Rastringin (PARTICLE position) {
	int i;
	long double fitness = 0.0;
	float pi = 3.1415;
	for (i = 0;   i < numberVariables; ++i)
		fitness = fitness+(position[i]*position[i]-10*cos(2.0*pi*position[i])+10);
	return fitness;
}

long double Griewank (PARTICLE position) {
	int i;
	long double fitness1, fitness2, fitness;
	fitness1 = 0.0;
	fitness2 = 1.0;
	 for (i = 1;    i < numberVariables+1; ++i) {
		 fitness1 = fitness1+(long double)(position[i]*position[i]);
		 fitness2 = fitness2*(cos(position[i]/sqrt((float)i)));
	 }
	 fitness = 1+(fitness1/4000)-fitness2;
	 return fitness;
}

long double Schaffer (PARTICLE position) {
	long double x, y;
	long double temp1, temp2;
	 x = (long double)position[0];
	 y = (long double)position[1];
	 temp1 = (long double)sin(sqrt(x * x + y * y));
	 temp2 = (long double)(1 + 0.001 * (x * x + y * y));
	 return (long double)(0.5 + (temp1 * temp1 - 0.5) / (temp2 * temp2));
}

long double weierstrass (PARTICLE position) {
	int i, j;
	long double res;
	long double sum;
	long double a, b;
	int k_max;
	a = 0.5;
	b = 3.0;
	k_max = 20;
	res = 0.0;
	for (i = 0; i < numberVariables; i++) {
		sum = 0.0;
		for (j = 0; j <= k_max; j++)
			sum = sum+pow(a,j)*cos(2.0*3.1415*pow(b,j)*(position[i]+0.5));
		res = res+sum;
	}
	sum = 0;
	//for (j = 0; j <= k_max; j++)
	//	sum = sum+pow(a,j)*cos(2.0*3.1415*pow(b,j)*0.5);
	res = res+60;
	return (res);
}

long double Schwefel (PARTICLE position) {
	int j;
	long double fitness;
	long double sum;
	sum = 0.0;
	for (j = 0; j < numberVariables; j++)
		sum = sum+position[j]*sin(sqrt(fabs((float)position[j])));
	fitness = 418.9829*numberVariables-sum;
	return (fitness);
}

long double evaluate (PARTICLE position) {
	long double fitness;
	if (problem == 1)	fitness = Sphere(position);
	if (problem == 2)   fitness = Rosenbrock(position);
	if (problem == 3)	fitness = Rastringin(position);
	if (problem == 4)   fitness = Griewank(position);
	if (problem == 5)   fitness = Schaffer(position);
	if (problem == 6)   fitness = weierstrass(position);
	if (problem == 7)   fitness = Ackley(position);
	if (problem == 8)   fitness = Schwefel(position);
	if (problem == 9)   fitness = Hyper(position);
	if (problem == 10)   fitness = Quartic(position);
	return fitness;
}
