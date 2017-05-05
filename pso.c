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

/// Internal function for hashing PRNG seed with thread ID
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

static const char * pso_validate_params(PSO_PARAMS params) {

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
PSO * pso_new(PSO_PARAMS params, pso_func func, unsigned int seed) {

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
	memmove(&pso->params, &params, sizeof(PSO_PARAMS));

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

	// Initialize best position vectors
	pso->best_position =
		(double *) malloc(pso->params.nvars * sizeof(double));
	pso->best_position_so_far =
		(double *) malloc(pso->params.nvars * sizeof(double));

	// Initialize cells and particles
	pso->particles =
		(PSO_PARTICLE *) malloc(pso->popSize * sizeof(PSO_PARTICLE));
	for (i = 0; i < pso->popSize; ++i) {
		pso->particles[i].informants_best_position_so_far =
			(double *) malloc(pso->params.nvars * sizeof(double));
		pso->particles[i].best_position_so_far =
			(double *) malloc(pso->params.nvars * sizeof(double));
		pso->particles[i].position =
			(double *) malloc(pso->params.nvars * sizeof(double));
		pso->particles[i].velocity =
			(double *) malloc(pso->params.nvars * sizeof(double));
	}

	pso->cell =
		(unsigned int **) malloc(pso->params.max_x * sizeof(unsigned int));

	z = 0;
	for (i = 0; i < pso->params.max_x; ++i){
		pso->cell[i] =
			(unsigned int *) malloc(pso->params.max_y * sizeof(unsigned int));
		for (j = 0; j < pso->params.max_y; ++j){
			// Set cell as occupied (particle is the id)
			pso->cell[i][j] = z;
			// Set particle default position
			pso->particles[z].x = i;
			pso->particles[z].y = j;
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
			pso->particles[i].position[j] =
				rds_luniform(&pso->prng_states[0], xmin, xmax);

			// Initialize velocity for current variable of current particle
			pso->particles[i].velocity[j] =
				rds_luniform(
					&pso->prng_states[0], -pso->params.Xmax, pso->params.Xmax) *
				(0.5 - rds_luniform(&pso->prng_states[0], 0, 1.0));

			// Set best position so far as current position
			pso->particles[i].best_position_so_far[j] =
				pso->particles[i].position[j];

			// Set best informat so far as myself
			pso->particles[i].informants_best_position_so_far[j] =
				pso->particles[i].position[j];
		}

		// Determine fitness for current particle
		pso->particles[i].fitness =
			pso->evaluate(pso->particles[i].position, pso->params.nvars);

		// Set my own fitness as best fitness so far
		pso->particles[i].best_fitness_so_far =
			pso->particles[i].fitness;

		// Set me as the best informant so far
		pso->particles[i].informants_best_fitness_so_far =
			pso->particles[i].fitness;
	}

	// Initialize remaining PSO variables to their defaults
	pso->best_so_far = pso->particles[0].fitness;
	pso->best_so_far_id = 0;
	pso->evaluations = 0;
	pso->minFitness = pso->particles[0].fitness;
	pso->worst_id = 0;

	// Return initialized PSO model
	return pso;

}

/**
 * Release a PSO model.
 *
 * @param[in] PSO model to release.
 */
void pso_destroy(PSO * pso) {

	// Release PRNG states
	free(pso->prng_states);

	// Release individual particles
	for (unsigned int i = 0; i < pso->popSize; ++i) {
		free(pso->particles[i].informants_best_position_so_far);
		free(pso->particles[i].best_position_so_far);
		free(pso->particles[i].position);
		free(pso->particles[i].velocity);
	}

	// Release particle vector
	free(pso->particles);

	// Free cells
	for (unsigned int i = 0; i < pso->params.max_y; ++i) {
		free(pso->cell[i]);
	}
	free(pso->cell);

	// Release best position vectors
	free(pso->best_position);
	free(pso->best_position_so_far);

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
	pso->best_fitness = pso->particles[0].fitness;
	pso->worst_fitness = 0;
	pso->average_fitness = 0;

	// Cycle through particles
	for (i = 0; i < pso->popSize; ++i) {

		// Updates worst in population
		if (pso->particles[i].fitness > pso->worst_fitness) {
			pso->worst_fitness = pso->particles[i].fitness;
			pso->worst_id = i;
		}

		// Updates best_so_far in population
		if (pso->particles[i].fitness < pso->best_so_far) {
			pso->best_so_far = pso->particles[i].fitness;
			memmove(pso->best_position_so_far,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));
		}

		// Updates best in current population
		if (pso->particles[i].fitness < pso->best_fitness) {
			pso->best_fitness = pso->particles[i].fitness;
			memmove(pso->best_position,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));
		}

		// Updates particle's best position
		if (pso->particles[i].fitness < pso->particles[i].best_fitness_so_far) {
			pso->particles[i].best_fitness_so_far = pso->particles[i].fitness;
			memmove(pso->particles[i].best_position_so_far,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));

		}

		// Updates best informant
		if (pso->particles[i].fitness <
				pso->particles[i].informants_best_fitness_so_far) {

			pso->particles[i].informants_best_fitness_so_far =
				pso->particles[i].fitness;

			memmove(pso->particles[i].informants_best_position_so_far,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));
		}

		pso->average_fitness = pso->average_fitness + pso->particles[i].fitness;
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
			pg = pso->particles[a].informants_best_position_so_far[i];
		if (pso->params.gbest == 1)
			// Global best
			pg = pso->best_position_so_far[i];

		pi = pso->particles[a].best_position_so_far[i];
		v = pso->particles[a].velocity[i];
		x = pso->particles[a].position[i];
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
		pso->particles[a].position[i] = x;
		pso->particles[a].velocity[i] = v;

	} // Cycle variables

	// Determine particle fitness for new position
	pso->particles[a].fitness =
		pso->evaluate(pso->particles[a].position, pso->params.nvars);

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

			// If a neighbor particle is the worst particle...
			if (neighParticle == pso->worst_id)
				// ...mark current particle for updating (SS-PSO only)
				update = 1;

			// Does the neighbor know of better fitness than current
			// particle?
			if (pso->particles[neighParticle].best_fitness_so_far <
					pso->particles[a].informants_best_fitness_so_far) {

				// If so, current particle will get that knowledge also
				pso->particles[a].informants_best_fitness_so_far =
					pso->particles[neighParticle].best_fitness_so_far;

				memmove(
					pso->particles[a].informants_best_position_so_far,
					pso->particles[neighParticle].best_position_so_far,
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
