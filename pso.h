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
extern const char * pso_error;

/// PSO model
typedef struct pso PSO;

/// Functions optimized by PSO
typedef double (* pso_func_opt)(double * vars, unsigned int nvars);

/// Hook functions
typedef void (* pso_func_hook)(PSO * pso);

// Structures

/// PSO parameters
typedef struct {
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

typedef struct {
	int dx;
	int dy;
} PSO_NEIGHBOR;

typedef struct {
	unsigned int num_neighs;
	const PSO_NEIGHBOR * neighs;
} PSO_NEIGHBORHOOD;

typedef struct {
	int x;
	int y;
	long double fitness;
	long double best_fitness_so_far;
	long double informants_best_fitness_so_far;
	double * informants_best_position_so_far;
	double * best_position_so_far;
	double * position;
	double * velocity;
} PSO_PARTICLE;

struct pso {
	PSO_PARAMS params;
	mt_state * prng_states;
	unsigned int popSize;
	unsigned int num_threads;
	pso_func_opt evaluate;

	unsigned int n_hooks;
	unsigned int alloc_hooks;
	pso_func_hook * hooks;

	/// Neighborhood used in PSO.
	const PSO_NEIGHBORHOOD * neighbors;

	double * best_position;
	double * best_position_so_far;

	PSO_PARTICLE * particles;
	unsigned int ** cell; // if ocupied, particle is the id, else, -1

	unsigned int iterations;
	unsigned int evaluations;
	unsigned int crit_evals;

	double minFitness;
	double average_fitness;
	double best_fitness;
	double best_so_far;
	double worst_fitness;
	//double worst_so_far;
	int best_so_far_id;
	int worst_id;
	//int worst_so_far_id;
};

// Initialize PSO model.
PSO * pso_new(PSO_PARAMS params, pso_func_opt func, unsigned int seed);

// Destroy PSO model.
void pso_destroy(PSO * pso);

// Update population data.
void pso_update_pop_data(PSO * pso);

// Update position and velocity of all or some of the particles.
void pso_update_particles(unsigned int iter, PSO * pso);

/// PSO algorithm
void pso_run(PSO * pso);

void pso_hook_add(PSO * pso, pso_func_hook func);


#endif
