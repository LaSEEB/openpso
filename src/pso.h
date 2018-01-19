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

#include "mtwist.h"

/// If PSO terminates with error, error message will be placed here
extern const char *pso_error;

/// PSO model
typedef struct pso PSO;

/// Forward declaration of PSO_PARTICLE
typedef struct pso_particle PSO_PARTICLE;

/// The PSO_NEIGH_ITERATOR can be anything, it is defined by the specific
/// topology instance
typedef void *PSO_NEIGH_ITERATOR;

/// Function which returns a neighbor iterator, defined by the specific
/// topology
typedef PSO_NEIGH_ITERATOR *(*pso_neigh_iterator)(PSO *, PSO_PARTICLE *);

/// Function which gets the next neighbor, defined by the specific topology
typedef PSO_PARTICLE *(*pso_neigh_next)(PSO_NEIGH_ITERATOR *);

/// Functions optimized by PSO
typedef double (*pso_func_opt)(double *vars, unsigned int nvars);

/// Hook functions
typedef void (*pso_func_hook)(PSO *pso);

// Structures

/// PSO parameters
typedef struct {

	// Initial number of particles
	unsigned int initPopSize;
	//
	unsigned int max_x, max_y;
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
	/// Stop criterion
	double crit;
	// Keep going until max_evaluations after stop criterion is meet?
	int crit_keep_going;

} PSO_PARAMS;

struct pso_particle {
	int x;
	int y;
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

	pso_neigh_iterator iterator;
	pso_neigh_next next;

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
