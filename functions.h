/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <malloc.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "definicoes.h"

typedef long double (*SelFunc)(PARTICLE position, unsigned int nvars);

long double Sphere(PARTICLE position, unsigned int nvars);
long double Rosenbrock(PARTICLE position, unsigned int nvars);
long double Rastrigin(PARTICLE position, unsigned int nvars);
long double Griewank(PARTICLE position, unsigned int nvars);
long double Ackley(PARTICLE position, unsigned int nvars);
long double Evaluate(PARTICLE position, unsigned int nvars);
long double weierstrass(PARTICLE position, unsigned int nvars);
long double Schwefel(PARTICLE position, unsigned int nvars);
long double Hyper(PARTICLE position, unsigned int nvars);
long double Quartic(PARTICLE position, unsigned int nvars);

SelFunc getSelFunc(unsigned int problem);
