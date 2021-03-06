/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Definitions for several Particle Swarm Optimization (PSO) algorithms.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

#ifndef __PSO_H_
#define __PSO_H_

#define PSO_PARAMS_MAXSTR 20

#include "common.h"
#include "mtwist.h"
#include "watershed.h"

/// If PSO terminates with error, error message will be placed here
extern const char *pso_error;

/// Forward declaration of PSO_PARTICLE
typedef struct pso_particle PSO_PARTICLE;

/// A PSO topology, can be anything
typedef void *PSO_TOPOLOGY;

typedef PSO_TOPOLOGY (*pso_topol_new)(PSO *);
typedef void (*pso_topol_destroy)(PSO_TOPOLOGY);

/// Function which restarts a neighbor iterator, defined by the specific
/// topology
typedef void (*pso_topol_iterate)(PSO_TOPOLOGY, PSO_PARTICLE *);

/// Function which gets the next neighbor, defined by the specific topology
typedef PSO_PARTICLE *(*pso_topol_next)(PSO_TOPOLOGY, PSO_PARTICLE *);

/// Functions optimized by PSO
typedef double (*pso_func_opt)(double *vars, unsigned int nvars);

/// Hook functions
typedef void (*pso_func_hook)(PSO *pso);

// Structures

/// PSO parameters
typedef struct {

	// Initial number of particles
	unsigned int initPopSize;

	// Topology parameters
	struct {
		pso_topol_new new;
		pso_topol_destroy destroy;
		pso_topol_iterate iterate;
		pso_topol_next next;
	} topol;

	//Maximum number of iterations
	unsigned int max_t;
	// Maximum number of evaluations
	unsigned int max_evaluations;
	// PSO algorithm to use
	unsigned int algorithm;
	//
	int gbest;
	//
	unsigned int neighborhood;
	//
	double Xmax;
	//
	double Vmax;
	//
	double chi;
	//
	double omega;
	//
	double c;
	//
	unsigned int nvars;
	//
	unsigned int iWeightStrategy, cStrategy;
	//
	int assyInitialization;
	//
	double initialXmin;
	//
	double initialXmax;
	// Number of extra random particles for small world PSO
	int numExtraRndNeighs;
	// Stop criterion
	double crit;
	// Keep going until max_evaluations after stop criterion is meet?
	int crit_keep_going;
	// Watershed strategy
	char watershed_strategy[PSO_PARAMS_MAXSTR];

} PSO_PARAMS;

struct pso_particle {
	//int x;
	//int y;
	void *neigh_info;

	unsigned int id;
	double fitness;
	double best_fitness_so_far;
	double informants_best_fitness_so_far;
	double *informants_best_position_so_far;
	double *best_position_so_far;
	double *position;
	double *velocity;
};

struct pso {

	PSO_PARAMS params;
	mt_state *prng_states;
	unsigned int pop_size;
	unsigned int num_threads;
	pso_func_opt evaluate;
	pso_watershed watershed;

	PSO_TOPOLOGY topol;

	unsigned int n_hooks;
	unsigned int alloc_hooks;
	pso_func_hook *hooks;

	double *best_position;
	double *best_position_so_far;

	PSO_PARTICLE *particles;

	unsigned int iterations;
	unsigned int evaluations;
	unsigned int crit_evals;

	double min_fitness;
	double average_fitness;
	double best_fitness;
	double best_so_far;
	double worst_fitness;
	double watershed_max_fitness;

	//double worst_so_far;
	PSO_PARTICLE *best_so_far_particle;
	PSO_PARTICLE *worst_particle;
	//int worst_so_far_id;
};

// Initialize PSO model.
PSO *pso_new(PSO_PARAMS params, pso_func_opt func, unsigned int seed);

// Destroy PSO model.
void pso_destroy(PSO *pso);

// Update population data.
void pso_update_pop_data(PSO *pso);

// Update position and velocity of all or some of the particles.
void pso_update_particles(unsigned int iter, PSO *pso);

// Run PSO algorithm
void pso_run(PSO *pso);

// Add end-of-iteration hook.
void pso_hook_add(PSO *pso, pso_func_hook func);

#endif
