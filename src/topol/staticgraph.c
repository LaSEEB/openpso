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
#include <string.h>

#define MAXFILELINE 256
#define MAXFILENAME 1024

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
	char tgf_file[MAXFILENAME];
	unsigned int nparticles;
	long file_edges_start;
} params;

unsigned int pso_staticgraph_parse_params(dictionary *d) {

	FILE *fp = NULL;
	char line[MAXFILELINE];
	int np = 0;
	strncpy(params.tgf_file,
		iniparser_getstring(d, "topology:tgf_file", "?"), MAXFILENAME);

	// Parse first part of TGF file, put configuration in params
	// Determine number of particles
	fp = fopen(params.tgf_file, "r");
	if (fp == NULL) {
		ERROR_EXIT("File not found: '%s'", params.tgf_file);
	}

	for (int i = 1; fgets(line, MAXFILELINE, fp) != NULL; ++i) {

		int nread;

		// Check if particle enumeration is finished
		if (line[0] == '#') break;

		// Check next particle
		nread = sscanf(line, "%d", &np);

		if (nread != 1) {
			ERROR_EXIT("Unexpected input at line %d of file '%s'",
				i, params.tgf_file);
		}

		if (np != i) {
			ERROR_EXIT("Expecting particle %d at line %d of file '%s' "
				"but got %d instead", i, i, params.tgf_file, np);
		}
	}

	if (line[0] != '#') {
		ERROR_EXIT("File '%s' ends too early, no line beginning with "
			"'#' found.", params.tgf_file);
	}

	params.file_edges_start = ftell(fp);
	fclose(fp);
	params.nparticles = np;


	return params.nparticles;
}

///////////// DURING PSO_NEW
PSO_TOPOLOGY pso_staticgraph_new(PSO *pso) {

	int np = params.nparticles;
	FILE *fp = NULL;
	char line[MAXFILELINE];
	PSO_SG_TOPOL *topol = malloc(sizeof(PSO_SG_TOPOL));

	topol->particles = pso->particles;

	// Initialize particles
	for (unsigned int i = 0; i < params.nparticles; ++i) {

		// Define neighbor information for this cell
		PSO_SG_NEIGH_INFO *ninfo = calloc(1, sizeof(PSO_SG_NEIGH_INFO));

		// Allocate memory for the first neighbor
		PSO_SG_ADJLIST *lst_elem = calloc(1, sizeof(PSO_SG_ADJLIST));

		// First neighbor is always the cell itself
		lst_elem->p = &topol->particles[i];
		lst_elem->next = NULL;
		ninfo->head = lst_elem;

		// Keep neighbor info for this cell
		topol->particles[i].neigh_info = ninfo;
	}

	// Parse second part of TGF file, determine edges between particles
	fp = fopen(params.tgf_file, "r");
	if (fp == NULL) {
		ERROR_EXIT("File not found: '%s'", params.tgf_file);
	}

	// Skip first part of file
	fseek(fp, params.file_edges_start, SEEK_SET);

	// Determine edges
	for (int i = np + 1; fgets(line, MAXFILELINE, fp) != NULL; ++i) {

		int nread;
		int p1, p2;
		PSO_SG_ADJLIST *lst_elem = NULL;
		PSO_SG_NEIGH_INFO *ninfo = NULL;

		// Check next particle
		nread = sscanf(line, "%d %d", &p1, &p2);

		if (nread != 2) {
			ERROR_EXIT("Unexpected input at line %d of file '%s'",
				i, params.tgf_file);
		}

		if ((p1 > np) || (p2 > np)) {
			ERROR_EXIT("Invalid particle ID at line %d of file '%s' ",
				i, params.tgf_file);
		}

		ninfo = (PSO_SG_NEIGH_INFO *) topol->particles[p1 - 1].neigh_info;
		lst_elem = calloc(1, sizeof(PSO_SG_ADJLIST));
		lst_elem->p = &topol->particles[p2 - 1];
		lst_elem->next = ninfo->head;
		ninfo->head = lst_elem;
	}

	return (PSO_TOPOLOGY) topol;
}


///////////// DURING PSO_DESTROY
void pso_staticgraph_destroy(PSO_TOPOLOGY topol) {

	PSO_SG_TOPOL *sgtopol = (PSO_SG_TOPOL *) topol;

	// Free neighbor info in each cell
	for (unsigned int i = 0; i < params.nparticles; ++i) {

		PSO_SG_NEIGH_INFO *ninfo =
			(PSO_SG_NEIGH_INFO *) sgtopol->particles[i].neigh_info;

		PSO_SG_ADJLIST *lst_elem = ninfo->head;

		// Free adjency list for each particle
		while (lst_elem != NULL) {
			PSO_SG_ADJLIST *aux = lst_elem->next;
			free(lst_elem);
			lst_elem = aux;
		}

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
