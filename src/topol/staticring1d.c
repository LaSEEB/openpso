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

#include "staticring1d.h"
#include "../errorhandling.h"
#include "zf_log.h"

typedef struct {
	unsigned int pos;
	unsigned int iter;
} PSO_SR1D_NEIGH_INFO;

typedef struct {
	unsigned int num_neighs;
	int *neighs;
	PSO_PARTICLE **particles;
} PSO_SR1D_TOPOL;

static struct {
	unsigned int nparticles;
	unsigned int radius;
} params;

unsigned int pso_staticring1d_parse_params(dictionary *d) {

	params.nparticles = (unsigned int)
			iniparser_getint(d, "topology:nparticles", 0);

	params.radius = (unsigned int)
			iniparser_getint(d, "topology:radius", 0);

	return params.nparticles;
}

///////////// DURING PSO_NEW
PSO_TOPOLOGY pso_staticring1d_new(PSO *pso) {

	PSO_SR1D_TOPOL *ring1d = malloc(sizeof(PSO_SR1D_TOPOL));

	// Setup neighbor mask

	// Initialize cells
	ring1d->particles = calloc(params.nparticles, sizeof(PSO_PARTICLE *));
	ring1d->num_neighs = 2 * params.radius;
	ring1d->neighs = malloc(ring1d->num_neighs * sizeof(int));

	for (unsigned int i = 0; i < ring1d->num_neighs; ++i) {
		int pos;
		if (i < ring1d->num_neighs / 2) {
			pos = i - (int) params.radius;
		} else {
			pos = i - (int) params.radius + 1;
		}
		ring1d->neighs[i] = pos;
	}

	for (unsigned int i = 0; i < params.nparticles; ++i) {

		// Define neighbor information for this cell
		PSO_SR1D_NEIGH_INFO *ninfo = malloc(sizeof(PSO_SR1D_NEIGH_INFO));
		ninfo->pos = i;
		ninfo->iter = 0;

		// Set cell as occupied
		ring1d->particles[i] = &pso->particles[i];

		// Keep neighbor info for this cell
		ring1d->particles[i]->neigh_info = ninfo;
	}

	return (PSO_TOPOLOGY) ring1d;
}


///////////// DURING PSO_DESTROY
void pso_staticring1d_destroy(PSO_TOPOLOGY topol) {

	PSO_SR1D_TOPOL *ring1d = (PSO_SR1D_TOPOL *) topol;

	free(ring1d->neighs);

	// Free cells
	for (unsigned int i = 0; i < params.nparticles; ++i) {
		free(ring1d->particles[i]->neigh_info);
	}
	free(ring1d->particles);
	free(ring1d);
}

/// Function which restarts a neighbor iterator, defined by the specific
/// topology
void pso_staticring1d_iterate(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	// Unused parameter
	(void)(topol);

	// Reset neighbor iterator for this particle
	PSO_SR1D_NEIGH_INFO *ninfo = (PSO_SR1D_NEIGH_INFO *) p->neigh_info;
	ninfo->iter = 0;
}

/// Function which gets the next neighbor, defined by the specific topology
PSO_PARTICLE *pso_staticring1d_next(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	PSO_SR1D_TOPOL *ring1d = (PSO_SR1D_TOPOL *) topol;
	PSO_SR1D_NEIGH_INFO *ninfo = (PSO_SR1D_NEIGH_INFO *) p->neigh_info;
	PSO_PARTICLE *neighParticle = NULL;

	if (ninfo->iter < ring1d->num_neighs) {

		int i = (int) ninfo->pos + ring1d->neighs[ninfo->iter];

		ninfo->iter++;

		// Adjust neighbors location according to ring topology
		if (i < 0)
			i = (int) params.nparticles - i;
		if (i >= (int) params.nparticles)
			i = i - (int) params.nparticles;

		// Get neighbor particle
		neighParticle = ring1d->particles[i];
	}

	return neighParticle;
}
