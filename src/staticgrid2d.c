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

#include "staticgrid2d.h"

// Known neighborhoods

/// Moore neighborhood
/*
static const PSO_NEIGHBORHOOD neighbors_moore = {
	.num_neighs = 9,
	.neighs = (PSO_NEIGHBOR[]) {{0, 0}, {0, 1}, {1, 1}, {1, 0},
		{1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}}};

/// Von Neumann neighborhood
static const PSO_NEIGHBORHOOD neighbors_vn = {
	.num_neighs = 5,
	.neighs = (PSO_NEIGHBOR[]) {{0, 0}, {0, 1}, {1, 0}, {0, -1}, {-1, 0}}};

/// Ring neighborhood (requires max_y == 1)
static const PSO_NEIGHBORHOOD neighbors_ring = {
	.num_neighs = 3,
	.neighs = (PSO_NEIGHBOR[]) {{-1, 0}, {0, 0}, {1, 0}}};
*/

/////////// PREVIOUSLY IN PSO.H, PSO object

unsigned int **cell; // if ocupied, particle is the id, else, -1

///////////// DURING PSO_NEW


unsigned int z;

// Setup neighbor mask
/*
if (pso->params.neighborhood == 0) { // Moore
	pso->neighbors = &neighbors_moore;
} else if (pso->params.neighborhood == 1) { // VN
	pso->neighbors = &neighbors_vn;
} else if (pso->params.neighborhood == 2) { // Ring (requires max_y == 1)
	pso->neighbors = &neighbors_ring;
	if (params.max_y != 1) {
		params.max_x = params.max_x * params.max_y;
		params.max_y = 1;
		ZF_LOGW("1D neighborhood selected, setting "
			"max_x==%u and max_y==%u\n", params.max_x, params.max_y);
	}
}


// Initialize cells
pso->cell =
	(unsigned int **) calloc(pso->params.max_x, sizeof(unsigned int *));

z = 0;
for (unsigned int i = 0; i < pso->params.max_x; ++i) {
	pso->cell[i] =
		(unsigned int *) calloc(pso->params.max_y, sizeof(unsigned int));
	for (unsigned int j = 0; j < pso->params.max_y; ++j) {
		// Set cell as occupied (particle is the id)
		pso->cell[i][j] = z;
		// Set particle default position
		pso->particles[z].x = i;
		pso->particles[z].y = j;
		// Increment id
		z++;
	}
}
*/

///////////// DURING PSO_DESTROY

/*
// Free cells
for (unsigned int i = 0; i < pso->params.max_y; ++i) {
	free(pso->cell[i]);
}
free(pso->cell);
*/

//////////// DURING PARTICLE UPDATE
/*
int i, j, ii, jj;
int neighParticle;

i = pso->particles[a].x + pso->neighbors->neighs[n].dx;
j = pso->particles[a].y + pso->neighbors->neighs[n].dy;

// Adjust neighbors location according to toroidal topology
ii = i;
jj = j;
if (i < 0)
	ii = pso->params.max_x - 1;
if (i >= (int) pso->params.max_x)
	ii = 0;
if (j < 0)
	jj = pso->params.max_y - 1;
if (j >= (int) pso->params.max_y)
	jj = 0;

// Get neighbor particle
neighParticle = pso->cell[ii][jj];
*/