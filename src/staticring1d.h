/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Topology stuff.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

#ifndef __PSO_STATICRING1D_H_
#define __PSO_STATICRING1D_H_

#include "pso.h"
#include "iniparser.h"

unsigned int pso_staticring1d_parse_params(dictionary *d);

PSO_TOPOLOGY pso_staticring1d_new(PSO *);
void pso_staticring1d_destroy(PSO_TOPOLOGY);

/// Function which restarts a neighbor iterator, defined by the specific
/// topology
void pso_staticring1d_iterate(PSO_TOPOLOGY, PSO_PARTICLE *);

/// Function which gets the next neighbor, defined by the specific topology
PSO_PARTICLE *pso_staticring1d_next(PSO_TOPOLOGY, PSO_PARTICLE *);


#endif
