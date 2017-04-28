/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "randistrs.h"
#include "mtwist.h"
#include "definicoes.h"
#include "functions_data.h"

typedef long double (* SelFunc)(PARTICLE position, unsigned int nvars, mt_state * prng_state);

long double Sphere(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Quadric(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Hyper(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Rastrigin(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Griewank(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Schaffer(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Weierstrass(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double Ackley(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double ShiftedQuadricWithNoise(PARTICLE position, unsigned int nvars, mt_state * prng_state);
long double RotatedGriewank(PARTICLE position, unsigned int nvars, mt_state * prng_state);

SelFunc getSelFunc(unsigned int problem);
