/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __PSO_H_
#define __PSO_H_

#include "mtwist.h"

#define MAX_POP_SIZE	8100
#define MAX_X			100
#define MAX_Y			100
#define MAX_VARIABLES   100

extern const char * pso_error;

typedef double (* pso_func)(double * vars, unsigned int nvars);


typedef struct {
	int dx;
	int dy;
} NEIGHBOR;

typedef struct {
	unsigned int num_neighs;
	const NEIGHBOR * neighs;
} NEIGHBORHOOD;

/////////////////////////////////////////////////////// STRUCTURES
typedef double PARTICLE[MAX_VARIABLES];

typedef struct {
	int x;
	int y;
	long double fitness;
	long double best_fitness_so_far;
	long double informants_best_fitness_so_far;
	PARTICLE informants_best_position_so_far;
	PARTICLE best_position_so_far;
	PARTICLE position;
	PARTICLE velocity;
} INDIVIDUAL;

typedef INDIVIDUAL SWARM[MAX_POP_SIZE];

typedef struct {
    int particle; // if ocupied, particle is the id, else, -1
} HABITAT;

typedef HABITAT CELL[MAX_X][MAX_Y];


// PSO parameters
typedef struct {
	//
	unsigned int max_x, max_y;
	//Maximum number of iterations
	unsigned int max_t;
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
} PARAMETERS;

typedef struct {
	PARAMETERS params;
	mt_state * prng_states;
	unsigned int popSize;
	unsigned int num_threads;
	pso_func evaluate;

	/// Neighborhood used in PSO.
	const NEIGHBORHOOD * neighbors;

    PARTICLE best_position;
	PARTICLE best_position_so_far;
	SWARM particle;
	CELL cell;
	long double minFitness;
	unsigned int evaluations;
	long double average_fitness;
	long double best_fitness;
    long double best_so_far;
	long double worst_fitness;
	//long double worst_so_far;
	int best_so_far_id;
    int worst_id;
	//int worst_so_far_id;
} PSO;

// Initialize PSO model.
PSO * pso_new(PARAMETERS params, pso_func func, unsigned int seed);

// Destroy PSO model.
void pso_destroy(PSO * pso);

// Update population data.
void updatePopulationData(PSO * pso);

// Update position and velocity of a particle.
void updateParticlePV(PSO * pso, int a, unsigned int iter);

// Update position and velocity of all or some of the particles.
void updateParticles(unsigned int iter, PSO * pso);


#endif
