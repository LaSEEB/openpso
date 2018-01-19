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

#ifndef __PSO_TOPOLOGY_H_
#define __PSO_TOPOLOGY_H_

#include "pso.h"

// TODO: Below stuff should go to specific grid2d topology plugin

typedef struct {
	int dx;
	int dy;
} PSO_NEIGHBOR;

typedef struct {
	unsigned int num_neighs;
	const PSO_NEIGHBOR *neighs;
} PSO_NEIGHBORHOOD;

#endif
