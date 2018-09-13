/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Definitions of several benchmark functions for optimization.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

#ifndef __FUNCTIONS_H_
#define __FUNCTIONS_H_

double Sphere(double *vars, unsigned int nvars);
double Quadric(double *vars, unsigned int nvars);
double Hyper(double *vars, unsigned int nvars);
double Rastrigin(double *vars, unsigned int nvars);
double Griewank(double *vars, unsigned int nvars);
double Schaffer6(double *vars, unsigned int nvars);
double Weierstrass(double *vars, unsigned int nvars);
double Ackley(double *vars, unsigned int nvars);
double ShiftedQuadricWithNoise(double *vars, unsigned int nvars);
double RotatedGriewank(double *vars, unsigned int nvars);
double Random01(double *vars, unsigned int nvars);

#endif
