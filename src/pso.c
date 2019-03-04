/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * An implementation of several Particle Swarm Optimization (PSO) algorithms.
 *
 * @author Carlos Fernandes
 * @author Nuno Fachada
 */

#ifdef _OPENMP
	#include <omp.h>
#endif
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <stdbool.h>

#include "pso.h"
#include "randistrs.h"
#include "zf_log.h"

// When hooks array is full/not initialized, by how much should we increment
// the hooks array capacity?
#define HOOKS_INC 2

/// If PSO terminates with error, error message will be placed here
const char *pso_error = NULL;

/// Variables with this structure will hold a fitness/particleID pair.
struct fit_id { double fit; unsigned int id; };

/**
 * Internal function for hashing PRNG seed with thread ID.
 *
 * @param[in] seed Global PRNG seed.
 * @param[in] tid Thread ID.
 * @return A thread-specific PRNG seed.
 */
static uint32_t pso_mixseed(uint32_t seed, uint32_t tid) {

	// A different base seed for each thread
	seed += tid;

	// Hash base seed
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

/**
 * Validate PSO parameters.
 *
 * @param[in] params The PSO parameters to validate.
 * @return `NULL` if no error is found, or the error string otherwise.
 */
static const char *pso_validate_params(PSO_PARAMS params) {

	// TODO Maybe initPopSize validation should be elsewhere?
	if (params.initPopSize < 1)
		return "Invalid initial population size";

	if (params.max_t < 1)
		return "Invalid input parameter: max_t";

	if (params.max_evaluations < 1)
		return "Invalid input parameter: max_evaluations";

	if ((params.algorithm < 1) || (params.algorithm > 2))
		return "Invalid input parameter: algorithm";

	if ((params.gbest < 0) || (params.gbest > 1))
		return "Invalid input parameter: gbest";

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

	if (params.iWeightStrategy > 1)
		return "Invalid input parameter: iWeightStrategy";

	if (params.cStrategy > 1)
		return "Invalid input parameter: cStrategy";

	if ((params.assyInitialization < 0) || (params.assyInitialization > 1))
		return "Invalid input parameter: assyInitialization";

	if (params.initialXmin < -DBL_MAX + 0.1) //FIXME
		return "Invalid input parameter: initialXmin";

	if (params.initialXmax < -DBL_MAX + 0.1) //FIXME
		return "Invalid input parameter: initialXmax";

	if (params.numExtraRndNeighs < 0)
		return "Invalid input parameter: numExtraRndNeighs";

	if (params.crit < -DBL_MAX + 0.1) //FIXME and am I needed here?
		return "Invalid input parameter: crit";

	if (params.crit_keep_going == -1) //FIXME and am I needed here?
		return "Invalid input parameter: crit_keep_going";

	return NULL;
}

/**
 * Update position and velocity of a particle.
 *
 * @param[in,out] PSO model containing particles to update.
 * @param[in] a Index of particle to update.
 * @param[in] iter Iteration count for current run.
 */
static void pso_update_particle_pv(PSO *pso, int a, unsigned int iter) {

	float v, x;
	float pi, pg = 0.0;
	double phi1, phi2;
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
	for (unsigned int i = 0; i < pso->params.nvars; ++i) {

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
 * Initialize PSO model.
 *
 * @param[in,out] PSO model to initialize.
 */
PSO *pso_new(PSO_PARAMS params, pso_func_opt func, unsigned int seed) {

	PSO *pso = NULL;
	float xmin, xmax;

	// Allocate memory for PSO
	pso = (PSO *) calloc(1, sizeof(PSO));

	// Validate parameters
	const char *error = pso_validate_params(params);
	if (error) {
		pso_error = error;
		ZF_LOGE("%s", error);
		return NULL;
	}

	// No hooks initially
	pso->n_hooks = 0;
	pso->alloc_hooks = 0;
	pso->hooks = NULL;

	// Keep the parameters and validate them
	memmove(&pso->params, &params, sizeof(PSO_PARAMS));

	// Determine the number of threads
	#ifdef _OPENMP
		pso->num_threads = omp_get_max_threads();
	#else
		pso->num_threads = 1;
	#endif

	// Initialize PRNG states
	pso->prng_states = (mt_state *) calloc(pso->num_threads, sizeof(mt_state));
	for (unsigned int i = 0; i < pso->num_threads; ++i) {
		mts_seed32new(&pso->prng_states[i], pso_mixseed(seed, i));
	}

	// Keep population size
	pso->pop_size = pso->params.initPopSize;

	// Keep function to evaluate
	pso->evaluate = func;

	// Initialize best position vectors
	pso->best_position_so_far =
		(double *) calloc(pso->params.nvars, sizeof(double));

	// Initialize particles vector
	pso->particles =
		(PSO_PARTICLE *) calloc(pso->pop_size, sizeof(PSO_PARTICLE));

	// Initialize individual particles
	for (unsigned int i = 0; i < pso->pop_size; ++i) {
		pso->particles[i].informants_best_position_so_far =
			(double *) calloc(pso->params.nvars, sizeof(double));
		pso->particles[i].best_position_so_far =
			(double *) calloc(pso->params.nvars, sizeof(double));
		pso->particles[i].position =
			(double *) calloc(pso->params.nvars, sizeof(double));
		pso->particles[i].velocity =
			(double *) calloc(pso->params.nvars, sizeof(double));
	}


	// Initialize topology ONLY AFTER particles are initialized
	pso->topol = pso->params.topol.new(pso);

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
	for (unsigned int i = 0; i < pso->pop_size; ++i) {

		// Initialize position and velocity for each variable of the current
		// particle
		for (unsigned int j = 0;  j < pso->params.nvars; ++j) {

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
	pso->best_so_far_particle = &pso->particles[0];
	pso->evaluations = 0;
	pso->min_fitness = pso->particles[0].fitness;
	pso->worst_particle = &pso->particles[0];

	// Return initialized PSO model
	return pso;

}

/**
 * Release a PSO model.
 *
 * @param[in] PSO model to release.
 */
void pso_destroy(PSO *pso) {

	// Release PRNG states
	free(pso->prng_states);

	// Release topology BEFORE releasing particles
	pso->params.topol.destroy(pso->topol);

	// Release individual particles
	for (unsigned int i = 0; i < pso->pop_size; ++i) {
		free(pso->particles[i].informants_best_position_so_far);
		free(pso->particles[i].best_position_so_far);
		free(pso->particles[i].position);
		free(pso->particles[i].velocity);
	}

	// Release particle vector
	free(pso->particles);

	// Release best position vectors
	free(pso->best_position_so_far);

	// Release hooks
	free(pso->hooks);

	// Release PSO model
	free(pso);
}

/**
 * Update population data. Let particles know about particle with best and
 * worst fitness in the population and calculate the average fitness in the
 * population.
 *
 * Client code will not usually call this function directly, unless more control
 * is desired in the way the PSO algorithm advances.
 *
 * @param[in,out] PSO model to update.
 */
void pso_update_pop_data(PSO *pso) {

	// Reset best and worst fitnesses
	struct fit_id
		worst_fitness = { .fit = 0.0, .id = 0 },
		best_fitness = { .fit = pso->particles[0].fitness, .id = 0 };
	double sum_fitness = 0;

// Uncomment the block below to parallelize the following for loop.
// Unfortunately, since this function executes really fast, the overhead of
// multithreading outweighs the potential parallelization gains.

/* #ifdef _OPENMP
	// Declare a "min" reduction operator which also keeps the particle ID
	#pragma omp declare reduction(min_fit_id : struct fit_id : \
		omp_out = omp_in.fit < omp_out.fit ? omp_in : omp_out) \
		initializer (omp_priv={ .fit = DBL_MAX, .id = 0 })
	// Declare a "max" reduction operator which also keeps the particle ID
	#pragma omp declare reduction(max_fit_id : struct fit_id : \
		omp_out = omp_in.fit > omp_out.fit ? omp_in : omp_out) \
		initializer (omp_priv={ .fit = 0, .id = 0 })
	// Parallel for pragma, perform reduction on worst_fitness, sum_fitness and
	// best_fitness
	#pragma omp parallel for \
		reduction(max_fit_id:worst_fitness) \
		reduction(+:sum_fitness) \
		reduction(min_fit_id:best_fitness)
#endif */

	// Cycle through particles
	for (unsigned int i = 0; i < pso->pop_size; ++i) {

		// Updates worst in population
		if (pso->particles[i].fitness > worst_fitness.fit) {

			worst_fitness.fit = pso->particles[i].fitness;
			worst_fitness.id = i;

		}

		// Updates best fitness/position in population (current iteration)
		if (pso->particles[i].fitness < best_fitness.fit) {

			best_fitness.fit = pso->particles[i].fitness;
			best_fitness.id = i;

		}

		// Updates current particle's best fitness/position so far (i.e. all
		// iterations so far)
		if (pso->particles[i].fitness < pso->particles[i].best_fitness_so_far) {

			pso->particles[i].best_fitness_so_far = pso->particles[i].fitness;

			memmove(pso->particles[i].best_position_so_far,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));

		}

		// Updates current particle's best informant so far (i.e. all iterations
		// so far)
		if (pso->particles[i].fitness <
				pso->particles[i].informants_best_fitness_so_far) {

			pso->particles[i].informants_best_fitness_so_far =
				pso->particles[i].fitness;

			memmove(pso->particles[i].informants_best_position_so_far,
				pso->particles[i].position,
				pso->params.nvars * sizeof(double));
		}

		sum_fitness += pso->particles[i].fitness;
	}

	// Update global worst fitness
	pso->worst_fitness = worst_fitness.fit;
	pso->worst_particle = &pso->particles[worst_fitness.id];

	// Update global best fitness/position for current iteration
	pso->best_fitness = best_fitness.fit;
	pso->best_position = pso->particles[best_fitness.id].position;

	// Updates best fitness/position so far in population (i.e. all
	// iterations so far)
	if (pso->best_so_far > best_fitness.fit) {

		pso->best_so_far = best_fitness.fit;
		pso->best_so_far_particle = &pso->particles[best_fitness.id];
		memmove(pso->best_position_so_far,
			pso->particles[best_fitness.id].position,
			pso->params.nvars * sizeof(double));
	}

	// Determine average fitness in the population
	pso->average_fitness = sum_fitness / pso->pop_size;
}

/**
 * Update position and velocity of all or some of the particles.
 *
 * Client code will not usually call this function directly, unless more control
 * is desired in the way the PSO algorithm advances.
 *
 * @param[in] iter Iteration count for current run.
 * @param[in,out] pso PSO model containing particles to move.
 */
void pso_update_particles(unsigned int iter, PSO *pso) {

	unsigned int evals = 0;

	// Cycle through particles
#ifdef _OPENMP
	unsigned int tid = omp_get_thread_num();
	#pragma omp parallel for reduction(+:evals) schedule(dynamic, 1)
#else
	unsigned int tid = 0;
#endif

	for (unsigned int a = 0; a < pso->pop_size; ++a) {

		// By default particle update is set to 0 (only relevant to SS-PSO)
		int update = 0;

		PSO_PARTICLE *currParticle = &pso->particles[a];
		PSO_PARTICLE *neighParticle = NULL;

		pso->params.topol.iterate(pso->topol, currParticle);

		// Cycle through neighbors
		while ((neighParticle =
			pso->params.topol.next(pso->topol, currParticle)) != NULL) {

			// If a neighbor particle is the worst particle...
			if (neighParticle == pso->worst_particle)
				// ...mark current particle for updating (SS-PSO only)
				update = 1;

			// Does the neighbor know of better fitness than current
			// particle?
			if (neighParticle->best_fitness_so_far <
					currParticle->informants_best_fitness_so_far) {

				// If so, current particle will get that knowledge also
				currParticle->informants_best_fitness_so_far =
					neighParticle->best_fitness_so_far;

				memmove(
					currParticle->informants_best_position_so_far,
					neighParticle->best_position_so_far,
					pso->params.nvars * sizeof(double));
			}

		} // Cycle neighbors

		// Consider extra neighbors (Small World PSO)?
		if (pso->params.numExtraRndNeighs > 0) {
			for (int i = 0; i < pso->params.numExtraRndNeighs; i++) {

				int rpID;
				PSO_PARTICLE* randomParticle;

				// Obtain random particle which is not a neighbor
				while (1) {

					// Obtain random particle
					rpID = rds_iuniform(
						&pso->prng_states[tid], 0, pso->pop_size);
					randomParticle = &pso->particles[rpID];

					// Is random particle the current particle?
					if (randomParticle == currParticle) {
						// If so, search for another random particle
						continue;
					}

					// Is random particle a neighbor?
					pso->params.topol.iterate(pso->topol, currParticle);
					while ((neighParticle =
						pso->params.topol.next(pso->topol, currParticle))
						!= NULL)
					{
						if (randomParticle == neighParticle)
						{
							// If so, search for another random particle
							continue;
						}
					}

					// If we get here, random particle is not a neighbor, so
					// break out of the loop
					break;
				}

				// Exchange knowledge with random particle

				// If a random particle is the worst particle...
				if (randomParticle == pso->worst_particle)
					// ...mark current particle for updating (SS-PSO only)
					update = 1;

				// Does the random particle know of better fitness than current
				// particle?
				if (randomParticle->best_fitness_so_far <
						currParticle->informants_best_fitness_so_far) {

					// If so, current particle will get that knowledge also
					currParticle->informants_best_fitness_so_far =
						randomParticle->best_fitness_so_far;

					memmove(
						currParticle->informants_best_position_so_far,
						randomParticle->best_position_so_far,
						pso->params.nvars * sizeof(double));
				}
			}
		}

		// Update current particle?
		if (
			// Regular PSO update strategy: always update
			(pso->params.algorithm == 1)
			// SS-PSO update strategy: only worst and its neighbors are updated
			|| ((pso->params.algorithm == 2) && (update == 1)))
		{
			// Update current particle
			pso_update_particle_pv(pso, a, iter);
			// Increment thread-local number of evaluations
			evals++;
		}

	 } // End parallel for - cycle particles

	 // Increment global number of evaluations
	 pso->evaluations += evals;

}

/**
 * Execute a complete PSO algorithm.
 *
 * Client code will usually call this function to run the PSO. For more control
 * client code can call the ::pso_update_pop_data() and ::pso_update_particles()
 * functions.
 *
 * @param[in,out] pso The PSO model object.
 */
void pso_run(PSO *pso) {

	// Aux. variables for current run
	pso->iterations = 0;
	pso->evaluations = 0;
	pso->crit_evals = 0;

	// Keep cycle going until maximum number of evaluations is reached
	do {

		// Update iteration count for current run
		pso->iterations++;

		// Let particles know about best and worst fitness and determine
		// average fitness
		pso_update_pop_data(pso);

		// Update all particles
		pso_update_particles(pso->iterations, pso);

		// Is the best so far below the stop criteria? If so did we already
		// saved the number of evaluations required to get below the
		// stop criteria?
		if ((pso->best_so_far < pso->params.crit) && (pso->crit_evals == 0)) {

			// Keep the number of evaluations which attained the stop criteria
			pso->crit_evals = pso->evaluations;

			// Stop current run if I'm not supposed to keep going
			if (!pso->params.crit_keep_going) break;

		}

		// Call end-of-iteration hook functions
		for (unsigned int i = 0; i < pso->n_hooks; ++i) {
			pso->hooks[i](pso);
		}

	} while (pso->evaluations < pso->params.max_evaluations);

}

/**
 * Add end-of-iteration hook function to the PSO model which will be call at
 * the end of each iteration of the complete PSO algorithm invoked with the
 * ::pso_run() function.
 *
 * @param[in,out] pso The PSO model object.
 * @param[in] func The hook function to be called at the end of each iteration.
 */
void pso_hook_add(PSO *pso, pso_func_hook func) {

	// Do we have space for more hooks?
	if (pso->n_hooks >= pso->alloc_hooks) {

		// If not, let's allocate some space
		ZF_LOGV("Allocating space for hooks. Was %u, now is %u",
			pso->alloc_hooks, pso->alloc_hooks + HOOKS_INC);

		pso->alloc_hooks = pso->alloc_hooks + HOOKS_INC;
		pso->hooks = (pso_func_hook *) realloc(pso->hooks,
			pso->alloc_hooks * sizeof(pso_func_hook));

	}

	// Add hook
	pso->hooks[pso->n_hooks] = func;
	pso->n_hooks++;

	ZF_LOGV("Registered hook, hooks count is %u/%u",
		pso->n_hooks, pso->alloc_hooks);

}
