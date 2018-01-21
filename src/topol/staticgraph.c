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

#include "staticgraph.h"
#include "../errorhandling.h"
#include "zf_log.h"

typedef struct pso_sg_adjlist PSO_SG_ADJLIST;

struct pso_sg_adjlist {
	PSO_PARTICLE *p;
	PSO_SG_ADJLIST *next;
};

typedef struct {
	PSO_SG_ADJLIST *head;
	PSO_SG_ADJLIST *current;
} PSO_SG_NEIGH_INFO;

typedef struct {
	PSO_PARTICLE *particles;
} PSO_SG_TOPOL;

static struct {
	const char *tgf_file;
	unsigned int nparticles;
} params;

unsigned int pso_staticgraph_parse_params(dictionary *d) {

	params.tgf_file = iniparser_getstring(d, "topology:tgf_file", 0);

	// Parse TGF file, put configuration in params
	// Determine number of particles

	return params.nparticles;
}

///////////// DURING PSO_NEW
PSO_TOPOLOGY pso_staticgraph_new(PSO *pso) {

	PSO_SG_TOPOL *topol = malloc(sizeof(PSO_SG_TOPOL));

	topol->particles = pso->particles;
	// Setup neighbor mask

	// Initialize cells
	for (unsigned int i = 0; i < params.nparticles; ++i) {

		// Define neighbor information for this cell
		PSO_SG_NEIGH_INFO *ninfo = malloc(sizeof(PSO_SG_NEIGH_INFO));

		// Read from whatever data structure we kept in params and create
		// adjacency lists for each particle

		// Keep neighbor info for this cell
		topol->particles[i].neigh_info = ninfo;
	}

	return (PSO_TOPOLOGY) topol;
}


///////////// DURING PSO_DESTROY
void pso_staticgraph_destroy(PSO_TOPOLOGY topol) {

	PSO_SG_TOPOL *sgtopol = (PSO_SG_TOPOL *) topol;

	// Free neighbor info in each cell
	for (unsigned int i = 0; i < params.nparticles; ++i) {

		// Need to free adjencey list for each particle

		free(sgtopol->particles[i].neigh_info);
	}

	free(topol);

}

/// Function which restarts a neighbor iterator, defined by the specific
/// topology
void pso_staticgraph_iterate(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	// Unused parameter
	(void)(topol);

	// Reset neighbor iterator for this particle
	PSO_SG_NEIGH_INFO *ninfo = (PSO_SG_NEIGH_INFO *) p->neigh_info;
	ninfo->current = ninfo->head;
}

/// Function which gets the next neighbor, defined by the specific topology
PSO_PARTICLE *pso_staticgraph_next(PSO_TOPOLOGY topol, PSO_PARTICLE *p) {

	(void)(topol);
	PSO_SG_NEIGH_INFO *ninfo = (PSO_SG_NEIGH_INFO *) p->neigh_info;
	PSO_PARTICLE *neighParticle = NULL;

	if (ninfo->current != NULL) {

		neighParticle = ninfo->current->p;
		ninfo->current = ninfo->current->next;
	}

	return neighParticle;
}
