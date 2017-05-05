/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef _OPENMP
	#include <omp.h>
#endif
#include <string.h>
#include <float.h>
#include <stdlib.h>

#include "pso.h"
#include "randistrs.h"
#include "zf_log.h"

const char * pso_error = NULL;

// Known neighborhoods

/// Moore neighborhood
static const NEIGHBORHOOD neighbors_moore = {
	.num_neighs = 9,
	.neighs = (NEIGHBOR[]) {{0, 0}, {0, 1}, {1, 1}, {1, 0},
		{1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}}};

/// Von Neumann neighborhood
static const NEIGHBORHOOD neighbors_vn = {
	.num_neighs = 5,
	.neighs = (NEIGHBOR[]) {{0, 0}, {0, 1}, {1, 0}, {0, -1}, {-1, 0}}};

/// Ring neighborhood (requires max_y == 1)
static const NEIGHBORHOOD neighbors_ring = {
	.num_neighs = 3,
	.neighs = (NEIGHBOR[]) {{-1, 0}, {0, 0}, {1, 0}}};

static uint32_t pso_mixseed(uint32_t seed, uint32_t tid) {

	// Hash seed
	seed  = (seed + 0x7ed55d16) + (seed << 12);
	seed  = (seed + 0x165667b1) + (seed << 5);
	seed  = (seed ^ 0xc761c23c) ^ (seed >> 19);
	seed  = (seed + 0xd3a2646c) ^ (seed << 9);
	seed  = (seed + 0xfd7046c5) + (seed << 3);
	seed  = (seed ^ 0xb55a4f09) ^ (seed >> 16);

	// Hash thread ID
	tid = (tid ^ 61) ^ (tid >> 16);
	tid = tid + (tid << 3);
	tid = tid ^ (tid >> 4);
	tid = tid * 0x27d4eb2d;
	tid = tid ^ (tid >> 15);

	// Return mixed seed
	return seed ^ tid;
}

static const char * pso_validate_params(PARAMETERS params) {

		if (params.max_x < 1)
			return "Invalid input parameter: max_x";

		if (params.max_y < 1)
			return "Invalid input parameter: max_y";

		if (params.max_t < 1)
			return "Invalid input parameter: max_t";

		if ((params.algorithm < 1) || (params.algorithm > 2))
			return "Invalid input parameter: algorithm";

		if ((params.gbest < 0) || (params.gbest > 1))
			return "Invalid input parameter: gbest";

		if (params.neighborhood > 2) { // Moore
			return "Invalid input parameter: neighborhood";
		}

		if (params.Xmax < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: Xmax";

		if (params.Vmax < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: Vmax";

		if (params.chi < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: chi";

		if (params.omega < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: omega";

		if (params.c < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: c";

		if (params.nvars < 1)
			return "Invalid input parameter: numberVariables";

		if (params.iWeightStrategy == 2) //FIXME
			return "Invalid input parameter: iWeightStrategy";

		if (params.cStrategy == 2) //FIXME
			return "Invalid input parameter: cStrategy";

		if (params.assyInitialization == -1) //FIXME
			return "Invalid input parameter: assyInitialization";

		if (params.initialXmin < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: initialXmin";

		if (params.initialXmax < -DBL_MAX + 0.1) //FIXME
			return "Invalid input parameter: initialXmax";

		if (params.crit < -DBL_MAX + 0.1) //FIXME and am I needed here?
			return "Invalid input parameter: crit";

		if (params.crit_keep_going == -1) //FIXME and am I needed here?
			return "Invalid input parameter: crit_keep_going";

		return NULL;
}

/**
 * Initialize PSO model.
 *
 * @param[in,out] PSO model to initialize.
 */
PSO * pso_new(PARAMETERS params, pso_func func, unsigned int seed) {

	PSO * pso;
	unsigned int i, j, z;
	float xmin, xmax;

	// Allocate memory for PSO
	pso = (PSO *) calloc(1, sizeof(PSO));

	// Validate parameters
	const char * error = pso_validate_params(params);
	if (error) {
		pso_error = error;
		ZF_LOGE("%s", error);
		return NULL;
	}

	// Keep the parameters and validate them
	memmove(&pso->params, &params, sizeof(PARAMETERS));

	// Setup neighbor mask
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

	// Determine the number of threads
	#ifdef _OPENMP
		pso->num_threads = omp_get_max_threads();
	#else
		pso->num_threads = 1;
	#endif

	// Initialize PRNG states
	pso->prng_states = (mt_state *) malloc(pso->num_threads * sizeof(mt_state));
	for (i = 0; i < pso->num_threads; ++i) {
		mts_seed32new(&pso->prng_states[i], pso_mixseed(seed, i));
	}

	// Keep population size
	pso->popSize = pso->params.max_x * pso->params.max_y;

	// Keep function to evaluate
	pso->evaluate = func;

	// Initialize cells and particles
	z = 0;
	for (i = 0; i < pso->params.max_x; ++i){
		for (j = 0; j < pso->params.max_y; ++j){
			// Set cell as occupied (particle is the id)
			pso->cell[i][j].particle = z;
			// Set particle default position
			pso->particle[z].x = i;
			pso->particle[z].y = j;
			// Increment id
			z = z + 1;
		}
	}

	// Set initial position bounds
	if (pso->params.assyInitialization == 1) {
		// Assymetric initialization of the population
		xmin = pso->params.initialXmin;
		xmax = pso->params.initialXmax;
	} else {
		// Normal initialization
		xmin = -pso->params.Xmax;
		xmax = pso->params.Xmax;
	}

	//Initialize position and velocity of each particle
	for (i = 0; i < pso->popSize; ++i) {

		// Initialize position and velocity for each variable of the current
		// particle
		for (j = 0;  j < pso->params.nvars; ++j) {

			// Initialize position for current variable of current particle
			pso->particle[i].position[j] =
				rds_luniform(&pso->prng_states[0], xmin, xmax);

			// Initialize velocity for current variable of current particle
			pso->particle[i].velocity[j] =
				rds_luniform(
					&pso->prng_states[0], -pso->params.Xmax, pso->params.Xmax) *
				(0.5 - rds_luniform(&pso->prng_states[0], 0, 1.0));

			// Set best position so far as current position
			pso->particle[i].best_position_so_far[j] =
				pso->particle[i].position[j];

			// Set best informat so far as myself
			pso->particle[i].informants_best_position_so_far[j] =
				pso->particle[i].position[j];
		}

		// Determine fitness for current particle
		pso->particle[i].fitness =
			pso->evaluate(pso->particle[i].position, pso->params.nvars);

		// Set my own fitness as best fitness so far
		pso->particle[i].best_fitness_so_far =
			pso->particle[i].fitness;

		// Set me as the best informant so far
		pso->particle[i].informants_best_fitness_so_far =
			pso->particle[i].fitness;
	}

	// Initialize remaining PSO variables to their defaults
	pso->best_so_far = pso->particle[0].fitness;
	pso->best_so_far_id = 0;
	pso->evaluations = 0;
	pso->minFitness = pso->particle[0].fitness;
	pso->worst_id = 0;

	return pso;

}

void pso_destroy(PSO * pso) {

	// Release PRNG states
	free(pso->prng_states);

	// Release PSO model
	free(pso);
}

/**
 * Update population data. Let particles know about particle with
 * best and worst fitness in the population and calculate the average fitness
 * in the population.
 *
 * @param[in,out] PSO model to update.
 */
void updatePopulationData(PSO * pso) {

	// Aux. variables
	unsigned int i;

	// Reset best and worst fitnesses
	pso->best_fitness = pso->particle[0].fitness;
	pso->worst_fitness = 0;
	pso->average_fitness = 0;

	// Cycle through particles
	for (i = 0; i < pso->popSize; ++i) {

		// Updates worst in population
		if (pso->particle[i].fitness > pso->worst_fitness) {
			pso->worst_fitness = pso->particle[i].fitness;
			pso->worst_id = i;
		}

		// Updates best_so_far in population
		if (pso->particle[i].fitness < pso->best_so_far) {
			pso->best_so_far = pso->particle[i].fitness;
			memmove(pso->best_position_so_far,
				pso->particle[i].position,
				pso->params.nvars * sizeof(double));
		}

		// Updates best in current population
		if (pso->particle[i].fitness < pso->best_fitness) {
			pso->best_fitness = pso->particle[i].fitness;
			memmove(pso->best_position,
				pso->particle[i].position,
				pso->params.nvars * sizeof(double));
		}

		// Updates particle's best position
		if (pso->particle[i].fitness < pso->particle[i].best_fitness_so_far) {
			pso->particle[i].best_fitness_so_far = pso->particle[i].fitness;
			memmove(pso->particle[i].best_position_so_far,
				pso->particle[i].position,
				pso->params.nvars * sizeof(double));

		}

		// Updates best informant
		if (pso->particle[i].fitness <
				pso->particle[i].informants_best_fitness_so_far) {

			pso->particle[i].informants_best_fitness_so_far =
				pso->particle[i].fitness;

			memmove(pso->particle[i].informants_best_position_so_far,
				pso->particle[i].position,
				pso->params.nvars * sizeof(double));
		}

		pso->average_fitness = pso->average_fitness + pso->particle[i].fitness;
	}

	// Determine average fitness in the population
	pso->average_fitness = pso->average_fitness / pso->popSize;
}

/**
 * Update position and velocity of a particle.
 *
 * @param[in,out] PSO model containing particles to update.
 * @param[in] a Index of particle to update.
 * @param[in] iter Iteration count for current run.
 */
void updateParticlePV(PSO * pso, int a, unsigned int iter) {

	unsigned int i;
	float v, x;
	float pi, pg = 0.0;
	long double phi1, phi2;
	float c1, c2;
	float maxIW, minIW;
	double omega;
#ifdef _OPENMP
	unsigned int tid = omp_get_thread_num();
#else
	unsigned int tid = 0;
#endif

	c1 = pso->params.c;
	c2 = pso->params.c;
	maxIW = 0.9;
	minIW = 0.4;


	// TVIW-PSO
	if (pso->params.iWeightStrategy == 1) {
		omega = minIW + (maxIW - minIW) *
			(((float) pso->params.max_t -
			(float) iter) / (float) pso->params.max_t);
	} else {
		omega = pso->params.omega;
	}

	//TVAC-PSO
	if (pso->params.cStrategy == 1) {
		c1 = (0.5 - 2.5) * ((float) iter / (float) pso->params.max_t) + 2.5;
		c2 = (2.5 - 0.5) * ((float) iter / (float) pso->params.max_t) + 0.5;
	}

	// Cycle through variables
	for (i = 0; i < pso->params.nvars; ++i) {

		// Use local best or global best?
		if (pso->params.gbest == 0)
			// Local best
			pg = pso->particle[a].informants_best_position_so_far[i];
		if (pso->params.gbest == 1)
			// Global best
			pg = pso->best_position_so_far[i];

		pi = pso->particle[a].best_position_so_far[i];
		v = pso->particle[a].velocity[i];
		x = pso->particle[a].position[i];
		phi1 = rds_luniform(&pso->prng_states[tid], 0.0, c1);
		phi2 = rds_luniform(&pso->prng_states[tid], 0.0, c2);

		// Determine updated Velocity
		v = omega * v + (float) phi1 * (pi - x) + (float) phi2 * (pg - x);
		if (v > pso->params.Vmax) v = pso->params.Vmax;
		if (v < -pso->params.Vmax) v = -pso->params.Vmax;

		// Determine updated Position
		x = x + v;
		if (x > pso->params.Xmax) {
			x = pso->params.Xmax;
			v = 0;
		}
		if (x < -pso->params.Xmax) {
			x = -pso->params.Xmax;
			v = 0;
		}

		// Update Position and Velocity for current variable
		pso->particle[a].position[i] = x;
		pso->particle[a].velocity[i] = v;

	} // Cycle variables

	// Determine particle fitness for new position
	pso->particle[a].fitness =
		pso->evaluate(pso->particle[a].position, pso->params.nvars);

}

/**
 * Update position and velocity of all or some of the particles.
 *
 * @param[in] iter Iteration count for current run.
 * @param[in,out] pso PSO model containing particles to move.
 */
void updateParticles(unsigned int iter, PSO * pso) {

	unsigned int evals = 0;

	// Cycle through particles
#ifdef _OPENMP
	#pragma omp parallel for reduction(+:evals)
#endif
	for (unsigned int a = 0; a < pso->popSize; ++a) {
		// By default particle update is set to 0 (only relevant to SS-PSO)
		int update = 0;

		// Cycle through neighbors
		for (unsigned int n = 0; n < pso->neighbors->num_neighs; ++n) {

			int i, j, ii, jj;
			int neighParticle;

			i = pso->particle[a].x + pso->neighbors->neighs[n].dx;
			j = pso->particle[a].y + pso->neighbors->neighs[n].dy;

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
			neighParticle = pso->cell[ii][jj].particle;

			// If a neighbor particle is the worst particle...
			if (neighParticle == pso->worst_id)
				// ...mark current particle for updating (SS-PSO only)
				update = 1;

			// Does the neighbor know of better fitness than current
			// particle?
			if (pso->particle[neighParticle].best_fitness_so_far <
					pso->particle[a].informants_best_fitness_so_far) {

				// If so, current particle will get that knowledge also
				pso->particle[a].informants_best_fitness_so_far =
					pso->particle[neighParticle].best_fitness_so_far;

				memmove(
					pso->particle[a].informants_best_position_so_far,
					pso->particle[neighParticle].best_position_so_far,
					pso->params.nvars * sizeof(double));
			}

		} // Cycle neighbors

		// Update current particle?
		if (
			// Regular PSO update strategy: always update
			(pso->params.algorithm == 1)
			// SS-PSO update strategy: only worst and its neighbors are updated
			|| ((pso->params.algorithm == 2) && (update == 1)))
		{
			// Update current particle
			updateParticlePV(pso, a, iter);
			// Increment thread-local number of evaluations
			evals++;
		}

	 } // End parallel for - cycle particles

	 // Increment global number of evaluations
	 pso->evaluations += evals;

}
