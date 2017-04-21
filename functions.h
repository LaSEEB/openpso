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

long double Sphere (PARTICLE position);
long double Rosenbrock (PARTICLE position);
long double Rastrigin (PARTICLE position);
long double Griewank (PARTICLE position);
long double Ackley (PARTICLE position);
long double Evaluate (PARTICLE position);
long double weierstrass (PARTICLE position);
long double Schwefel (PARTICLE position);
long double Hyper (PARTICLE position);
long double Quartic (PARTICLE position);
long double evaluate (PARTICLE position);
