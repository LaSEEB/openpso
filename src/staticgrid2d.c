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
#include "iniparser.h"
#include "errorhandling.h"
#include "zf_log.h"

typedef struct {
	int dx;
	int dy;
} PSO_NEIGHBOR;

typedef struct {
	unsigned int num_neighs;
	const PSO_NEIGHBOR *neighs;
} PSO_NEIGHBORHOOD;

typedef struct {
	const PSO_NEIGHBORHOOD *nhood;
	unsigned int **cell; // if ocupied, particle is the id, else, -1
} PSO_GRID;

static struct {
	enum { MOORE, VON_NEUMANN, RING } neighborhood;
	unsigned int xdim;
	unsigned int ydim;
} params;

// Known neighborhoods

/// Moore neighborhood
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


unsigned int pso_staticgrid2d_parse_params(dictionary *d) {

	const char *neighborhood =
		iniparser_getstring(d, "topology:neighborhood", "Moore");

	if (strcmp(neighborhood, "MOORE") == 0) {
		params.neighborhood = MOORE;
	} else if (strcmp(neighborhood, "VN") == 0) {
		params.neighborhood = VON_NEUMANN;
	} else {
		ERROR_EXIT("Unknown Grid2D neighborhood: '%s'", neighborhood);
	}

	params.xdim = (unsigned int)
			iniparser_getint(d, "topology:xdim", 0);

	params.ydim = (unsigned int)
			iniparser_getint(d, "topology:ydim", 0);

	return params.xdim * params.ydim;
}

///////////// DURING PSO_NEW
void *pso_staticgrid2d_new() {

    PSO_GRID *grid2d = malloc(sizeof(PSO_GRID));

	unsigned int z = 0;

	// Setup neighbor mask

	if (params.neighborhood == MOORE) { // Moore
		grid2d->nhood = &neighbors_moore;
	} else if (params.neighborhood == VON_NEUMANN) { // VN
		grid2d->nhood = &neighbors_vn;
	} else if (params.neighborhood == RING) { // Ring (requires max_y == 1)
		grid2d->nhood = &neighbors_ring;
		if (params.ydim != 1) {
			params.xdim = params.xdim * params.ydim;
			params.ydim = 1;
			ZF_LOGW("1D neighborhood selected, setting "
				"xdim==%u and ydim==%u\n", params.xdim, params.ydim);
		}
	}

	// Initialize cells
	grid2d->cell =
		(unsigned int **) calloc(params.xdim, sizeof(unsigned int *));

	for (unsigned int i = 0; i < params.xdim; ++i) {
		grid2d->cell[i] =
			(unsigned int *) calloc(params.ydim, sizeof(unsigned int));
		for (unsigned int j = 0; j < params.ydim; ++j) {
			// Set cell as occupied (particle is the id)
			grid2d->cell[i][j] = z;
			// Increment id
			z++;
		}
	}

	return (void *) grid2d;
}


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
