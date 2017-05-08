/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Definition of a no-seed uniform PRNG function which returns doubles between
 * 0 and 1. The PRNG is used to create noise for some of the benchmarking
 * optimization functions.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

#ifndef __FPRNG_H_
#define __FPRNG_H_

// Return random double between 0 and 1.
double runif01(double * vars, unsigned int nvars);

#endif
