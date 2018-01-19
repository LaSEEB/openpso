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
	unsigned int xpos;
	unsigned int ypos;
	unsigned int iter;
} PSO_SG2D_NEIGH_INFO;

typedef struct {
	const PSO_NEIGHBORHOOD *nhood;
	PSO_PARTICLE ***particles; // if not ocupied position set to NULL
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
PSO_TOPOLOGY pso_staticgrid2d_new(PSO *pso) {

	PSO_GRID *grid2d = malloc(sizeof(PSO_GRID));

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
	grid2d->particles = calloc(params.xdim, sizeof(PSO_PARTICLE **));

	for (unsigned int x = 0; x < params.xdim; ++x) {
		grid2d->particles[x] = calloc(params.ydim, sizeof(PSO_PARTICLE *));
		for (unsigned int y = 0; y < params.ydim; ++y) {

			// Define neighbor information for this cell
			PSO_SG2D_NEIGH_INFO *ninfo = malloc(sizeof(PSO_SG2D_NEIGH_INFO));
			ninfo->xpos = x;
			ninfo->ypos = y;
			ninfo->iter = 0;

			// Set cell as occupied
			grid2d->particles[x][y] = &pso->particles[y * params.xdim + x];

			// Keep neighbor info for this cell
			grid2d->particles[x][y]->neigh_info	= ninfo;
		}
	}

	// Setup neighbor iterator functions
	pso->iterate = pso_grid2d_iterate;
	pso->next = pso_grid2d_next;

	return (void *) grid2d;
}


///////////// DURING PSO_DESTROY
void pso_staticgrid2d_destroy(PSO_TOPOLOGY topol) {

	PSO_GRID *grid2d = (PSO_GRID *) topol;

	// Free cells
	for (unsigned int x = 0; x < params.xdim; ++x) {
		for (unsigned int y = 0; y < params.ydim; ++y) {
			free(grid2d->particles[x][y]->neigh_info);
		}
		free(grid2d->particles[x]);
	}
	free(grid2d->particles);
}

/// Function which restarts a neighbor iterator, defined by the specific
/// topology
void pso_grid2d_iterate(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	// Unused parameter
	(void)(topol);

	// Reset neighbor iterator for this particle
	PSO_SG2D_NEIGH_INFO *ninfo = (PSO_SG2D_NEIGH_INFO *) p->neigh_info;
	ninfo->iter = 0;
}

/// Function which gets the next neighbor, defined by the specific topology
PSO_PARTICLE *pso_grid2d_next(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	PSO_GRID *grid2d = (PSO_GRID *) topol;
	PSO_SG2D_NEIGH_INFO *ninfo = (PSO_SG2D_NEIGH_INFO *) p->neigh_info;
	PSO_PARTICLE *neighParticle = NULL;

	int i, j, ii, jj;

	if (ninfo->iter < grid2d->nhood->num_neighs) {
		i = ninfo->xpos + grid2d->nhood->neighs[ninfo->iter].dx;
		j = ninfo->ypos + grid2d->nhood->neighs[ninfo->iter].dy;

		ninfo->iter++;

		// Adjust neighbors location according to toroidal topology
		ii = i;
		jj = j;
		if (i < 0)
			ii = params.xdim - 1;
		if (i >= (int) params.xdim)
			ii = 0;
		if (j < 0)
			jj = params.ydim - 1;
		if (j >= (int) params.ydim)
			jj = 0;

		// Get neighbor particle
		neighParticle = grid2d->particles[ii][jj];
	}

	return neighParticle;
}
