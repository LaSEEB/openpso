/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 /**
  * @file
  * Implementation of Steady-State Particle Swarm Optimization (SSPSO).
  *
  * @author Carlos Fernandes
  * @author Nuno Fachada
  */

// System libraries
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

// Local libraries
#include "randistrs.h"
#include "mtwist.h"
#include "iniparser.h"

// SSPSO headers
#include "definicoes.h"
#include "functions.h"

// Constants
#define DEFAULT_INPUT_FILE "input.ini"
#define DEFAULT_PRNG_SEED 1234
#define MAX_FILENAME_LEN 255

// *************** Global Variables *****************

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
	.neighs = (NEIGHBOR[]) {{-1, 0}, {0, 0}, {0, 1}}};

/// Seed for pseudo-random number generator.
static unsigned int prng_seed;
/// File containing PSO parameters.
static char input_file[MAX_FILENAME_LEN];
/// Function/problem to solve.
static SelFunc evaluate;
/// Neighborhood used in PSO.
static const NEIGHBORHOOD * neighbors;
/// Population size.
static unsigned int popSize;
/// Number of threads.
static unsigned int num_threads;
/// PRNG states
static mt_state * prng_states;

// PSO parameters

//
static unsigned int max_x, max_y;
//Maximum number of iterations
static unsigned int max_t;
// Maximum number of evaluations
static unsigned int max_evaluations;
// Number of runs
static unsigned int n_runs;
// PSO algorithm to use
static unsigned int algorithm;
//
static int gbest;
//
static unsigned int neighborhood;
//
static unsigned int problem;
//
static double Xmax;
//
static double Vmax;
//
static double chi;
//
static double omega;
//
static double c;
//
static unsigned int numberVariables;
//
static unsigned int iWeightStrategy, cStrategy;
//
static int assyInitialization;
//
static double initialXmin;
//
static double initialXmax;
/// Stop criterion
static double crit;
// Keep going until max_evaluations after stop criterion is meet?
static int crit_keep_going;

/**
 * Parse command-line options and read PSO parameters from file.
 *
 * @param[in] argc Number of program arguments.
 * @param[in] argv Program arguments. Argument at index 1 is the file containing
 * the PSO parameters. Argument at index 2 is the PRNG seed.
 * @return `EXIT_SUCCESS` if program executes successfully, `EXIT_FAILURE`
 * otherwise.
 */
static void parse_params(int argc, char * argv[]) {

	// INI object
	dictionary * ini;

	// Did user specify a PSO parameter file?
	if (argc >= 2)
		strncpy(input_file, argv[1], MAX_FILENAME_LEN);
	else
		strncpy(input_file, DEFAULT_INPUT_FILE, MAX_FILENAME_LEN);

	// Did user specify a PRNG seed?
	if (argc >= 3)
		prng_seed = atoi(argv[2]);
	else
		prng_seed = DEFAULT_PRNG_SEED;

	// Number of threads
#ifdef _OPENMP
	num_threads = omp_get_max_threads();
#else
	num_threads = 1;
#endif
	printf("Num threads=%u\n", num_threads);

	// Try to open PSO parameters file
	ini = iniparser_load(input_file);
	if (!ini) ERROR_EXIT("Unable to parse input file '%s'\n", input_file);

	// Read PSO parameters file

	max_x = (unsigned int) iniparser_getint(ini, "pso:max_x", 0);
	if (max_x < 1)
		ERROR_EXIT("Invalid input parameter: max_x\n");

	max_y = (unsigned int) iniparser_getint(ini, "pso:max_y", 0);
	if (max_y < 1)
		ERROR_EXIT("Invalid input parameter: max_y\n");

	n_runs = (unsigned int) iniparser_getint(ini, "pso:n_runs", 0);
	if (n_runs < 1)
		ERROR_EXIT("Invalid input parameter: n_runs\n");

	max_t = (unsigned int) iniparser_getint(ini, "pso:max_t", 0);
	if (max_t < 1)
		ERROR_EXIT("Invalid input parameter: max_t\n");

	max_evaluations = (unsigned int) iniparser_getint(ini, "pso:max_evaluations", 0);
	if (max_evaluations < 1)
		ERROR_EXIT("Invalid input parameter: max_evaluations\n");

	algorithm = (unsigned int) iniparser_getint(ini, "pso:algorithm", 0);
	if ((algorithm < 1) || (algorithm > 2))
		ERROR_EXIT("Invalid input parameter: algorithm\n");

	gbest = iniparser_getboolean(ini, "pso:gbest", -1);
	if (gbest == -1)
		ERROR_EXIT("Invalid input parameter: gbest\n");

	// Neighborhood
	neighborhood = (unsigned int) iniparser_getint(ini, "pso:neighborhood", 3);
	if (neighborhood == 0) { // Moore
		neighbors = &neighbors_moore;
	} else if (neighborhood == 1) { // VN
		neighbors = &neighbors_vn;
	} else if (neighborhood == 2) { // Ring (requires max_y == 1)
		neighbors = &neighbors_ring;
		if (max_y != 1) ERROR_EXIT("Ring neighborhood requires max_y==1\n");
	} else {
		ERROR_EXIT("Invalid input parameter: neighborhood\n");
	}

	problem = (unsigned int) iniparser_getint(ini, "pso:problem", 0);
	if ((problem < 1) || (problem > 10))
		ERROR_EXIT("Invalid input parameter: problem\n");

	Xmax = iniparser_getdouble(ini, "pso:xmax", -DBL_MAX);
	if (Xmax < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: Xmax\n");

	Vmax = iniparser_getdouble(ini, "pso:vmax", -DBL_MAX);
	if (Vmax < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: Vmax\n");

	chi = iniparser_getdouble(ini, "pso:chi", -DBL_MAX);
	if (chi < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: chi\n");

	omega = iniparser_getdouble(ini, "pso:omega", -DBL_MAX);
	if (omega < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: omega\n");

	c = iniparser_getdouble(ini, "pso:c", -DBL_MAX);
	if (c < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: c\n");

	numberVariables = (unsigned int) iniparser_getint(ini, "pso:numbervariables", 0);
	if (numberVariables < 1)
		ERROR_EXIT("Invalid input parameter: numberVariables\n");

	iWeightStrategy = (unsigned int) iniparser_getint(ini, "pso:iweightstrategy", 2);
	if (iWeightStrategy == 2)
		ERROR_EXIT("Invalid input parameter: iWeightStrategy\n");

	cStrategy = (unsigned int) iniparser_getint(ini, "pso:cstrategy", 2);
	if (cStrategy == 2)
		ERROR_EXIT("Invalid input parameter: cStrategy\n");

	assyInitialization = iniparser_getboolean(ini, "pso:assyinitialization", -1);
	if (assyInitialization == -1)
		ERROR_EXIT("Invalid input parameter: assyInitialization\n");

	initialXmin = iniparser_getdouble(ini, "pso:initialxmin", -DBL_MAX);
	if (initialXmin < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: initialXmin\n");

	initialXmax = iniparser_getdouble(ini, "pso:initialxmax", -DBL_MAX);
	if (initialXmax < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: initialXmax\n");

	crit = iniparser_getdouble(ini, "pso:crit", -DBL_MAX);
	if (crit < -DBL_MAX + 0.1)
		ERROR_EXIT("Invalid input parameter: crit\n");

 	crit_keep_going = iniparser_getboolean(ini, "pso:crit_keep_going", -1);
	if (crit_keep_going == -1)
		ERROR_EXIT("Invalid input parameter: crit_keep_going\n");

	// Determine population size based on grid size
	popSize = max_x * max_y;

	// Release dictionary object
	iniparser_freedict(ini);

}

/**
 * Initialize PSO model.
 *
 * @param[in,out] PSO model to initialize.
 */
static void initialize(MODEL * pso) {

	unsigned int i, j, z;
	float xmin, xmax;

	// Initialize cells and particles
	z = 0;
	for (i = 0; i < max_x; ++i){
		for (j = 0; j < max_y; ++j){
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
	if (assyInitialization == 1) {
		// Assymetric initialization of the population
		xmin = initialXmin;
		xmax = initialXmax;
	} else {
		// Normal initialization
		xmin = -Xmax;
		xmax = Xmax;
	}

	//Initialize position and velocity of each particle
	for (i = 0; i < popSize; ++i) {

		// Initialize position and velocity for each variable of the current
		// particle
		for (j = 0;  j < numberVariables; ++j) {

			// Initialize position for current variable of current particle
			pso->particle[i].position[j] =
				rds_luniform(&prng_states[0], xmin, xmax);

			// Initialize velocity for current variable of current particle
			pso->particle[i].velocity[j] =
				rds_luniform(&prng_states[0], -Xmax, Xmax) *
				(0.5 - rds_luniform(&prng_states[0], 0, 1.0));

			// Set best position so far as current position
			pso->particle[i].best_position_so_far[j] =
				pso->particle[i].position[j];

			// Set best informat so far as myself
			pso->particle[i].informants_best_position_so_far[j] =
				pso->particle[i].position[j];
		}

		// Determine fitness for current particle
		pso->particle[i].fitness =
			evaluate(pso->particle[i].position, numberVariables, &prng_states[0]);

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

}

/**
 * Update population data. Let particles know about particle with
 * best and worst fitness in the population and calculate the average fitness
 * in the population.
 *
 * @param[in,out] PSO model to update.
 */
static void updatePopulationData(MODEL * pso) {

	// Aux. variables
	unsigned int i;

	// Reset best and worst fitnesses
	pso->best_fitness = pso->particle[0].fitness;
	pso->worst_fitness = 0;
	pso->average_fitness = 0;

	// Cycle through particles
	for (i = 0; i < popSize; ++i) {

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
				numberVariables * sizeof(float));
		}

		// Updates best in current population
		if (pso->particle[i].fitness < pso->best_fitness) {
			pso->best_fitness = pso->particle[i].fitness;
			memmove(pso->best_position,
				pso->particle[i].position,
				numberVariables * sizeof(float));
		}

		// Updates particle's best position
		if (pso->particle[i].fitness < pso->particle[i].best_fitness_so_far) {
			pso->particle[i].best_fitness_so_far = pso->particle[i].fitness;
			memmove(pso->particle[i].best_position_so_far,
				pso->particle[i].position,
				numberVariables * sizeof(float));

		}

		// Updates best informant
		if (pso->particle[i].fitness <
				pso->particle[i].informants_best_fitness_so_far) {

			pso->particle[i].informants_best_fitness_so_far =
				pso->particle[i].fitness;

			memmove(pso->particle[i].informants_best_position_so_far,
				pso->particle[i].position,
				numberVariables * sizeof(float));
		}

		pso->average_fitness = pso->average_fitness + pso->particle[i].fitness;
	}

	// Determine average fitness in the population
	pso->average_fitness = pso->average_fitness / popSize;
}

/**
 * Update position and velocity of a particle.
 *
 * @param[in,out] PSO model containing particles to update.
 * @param[in] a Index of particle to update.
 * @param[in] iter Iteration count for current run.
 */
static void updateParticlePV(MODEL * pso, int a, unsigned int iter) {

	unsigned int i;
	float v, x;
	float pi, pg = 0.0;
	long double phi1, phi2;
	float c1, c2;
	float maxIW, minIW;
#ifdef _OPENMP
	unsigned int tid = omp_get_thread_num();
#else
	unsigned int tid = 0;
#endif

	c1 = c;
	c2 = c;
	maxIW = 0.9;
	minIW = 0.4;


	// TVIW-PSO
	if (iWeightStrategy == 1) {
		omega = minIW + (maxIW - minIW) *
			(((float) max_t - (float) iter) / (float) max_t);
	}

	//TVAC-PSO
	if (cStrategy == 1) {
		c1 = (0.5 - 2.5) * ((float) iter / (float) max_t) + 2.5;
		c2 = (2.5 - 0.5) * ((float) iter / (float) max_t) + 0.5;
	}

	// Cycle through variables
	for (i = 0; i < numberVariables; ++i) {

		// Use local best or global best?
		if (gbest == 0)
			// Local best
			pg = pso->particle[a].informants_best_position_so_far[i];
		if (gbest == 1)
			// Global best
			pg = pso->best_position_so_far[i];

		pi = pso->particle[a].best_position_so_far[i];
		v = pso->particle[a].velocity[i];
		x = pso->particle[a].position[i];
		phi1 = rds_luniform(&prng_states[tid], 0.0, c1);
		phi2 = rds_luniform(&prng_states[tid], 0.0, c2);

		// Determine updated Velocity
		v = omega * v + (float) phi1 * (pi - x) + (float) phi2 * (pg - x);
		if (v > Vmax) v = Vmax;
		if (v < -Vmax) v = -Vmax;

		// Determine updated Position
		x = x + v;
		if (x > Xmax) {
			x = Xmax;
			v = 0;
		}
		if (x < -Xmax) {
			x = -Xmax;
			v = 0;
		}

		// Update Position and Velocity for current variable
		pso->particle[a].position[i] = x;
		pso->particle[a].velocity[i] = v;

	} // Cycle variables

	// Determine particle fitness for new position
	pso->particle[a].fitness =
		evaluate(pso->particle[a].position, numberVariables, &prng_states[tid]);

	// Increment number of evaluations
	pso->evaluations += 1;

}

/**
 * Update position and velocity of all or some of the particles.
 *
 * @param[in] iter Iteration count for current run.
 * @param[in,out] pso PSO model containing particles to move.
 */
static void updateParticles(unsigned int iter, MODEL * pso) {

	unsigned int a, n;
	int i, j, ii, jj;
	int neighParticle;
	int update;

	// Cycle through particles
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (a = 0; a < popSize; ++a) {

		// By default particle update is set to 0 (only relevant to SS-PSO)
		update = 0;

		// Cycle through neighbors
		for (n = 0; n < neighbors->num_neighs; ++n) {

			i = pso->particle[a].x + neighbors->neighs[n].dx;
			j = pso->particle[a].y + neighbors->neighs[n].dy;

			// Adjust neighbors location according to toroidal topology
			ii = i;
			jj = j;
			if (i < 0)
				ii = max_x - 1;
			if (i >= (int) max_x)
				ii = 0;
			if (j < 0)
				jj = max_y - 1;
			if (j >= (int) max_y)
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
					numberVariables * sizeof(float));
			}

		} // Cycle neighbors

		// Typical PSO update strategy: update all
		if (algorithm == 1)
			updateParticlePV(pso, a, iter);

		// SS-PSO update strategy: only the worst and its neighbors are updated
		if ((algorithm == 2) && (update == 1))
			updateParticlePV(pso, a, iter);

	 } // Cycle particles

}

/**
 * Program starts here.
 *
 * @param[in] argc Number of program arguments.
 * @param[in] argv Program arguments. Argument at index 1 is the file containing
 * the PSO parameters. Argument at index 2 is the PRNG seed.
 * @return `EXIT_SUCCESS` if program executes successfully, `EXIT_FAILURE`
 * otherwise.
 */
int main(int argc, char* argv[]) {

	// Parse command-line arguments and read PSO parameter file.
	parse_params(argc, argv);

	FILE * out;
	unsigned int i;
	MODEL * pso;

	// Iteration count for current run
	unsigned int iter;

	unsigned int counter = 0;
	long double averageBestSoFar[max_evaluations / 100];
	unsigned int crit_evals[n_runs];
	double best_so_far[n_runs];
	int flag;

	// Show parameter information to user
	printf("PSO parameter file : %s\n", input_file);
	printf("PRNG seed          : %d\n", prng_seed);

	// Initialize PRNG states
	prng_states = (mt_state *) malloc(num_threads * sizeof(mt_state));
	for (i = 0; i < num_threads; ++i)
		mts_seed32new(&prng_states[i], prng_seed + i);

	// Set function/problem to optimize
	evaluate = getSelFunc(problem);

	// Initialize averageBestSoFar array, set contents to zero
	memset(averageBestSoFar, 0, max_evaluations / 100 * sizeof(long double));

	// Perform PSO runs
	for (i = 0; i < n_runs; ++i) {

		// Create PSO model for current run, set contents to zero
		pso = (MODEL *) calloc(1, sizeof(MODEL));

		// Initialize PSO for current run
		initialize(pso);

		// Aux. variables for current run
		flag = 0;
		iter = 0;
		counter = 0;

		// PSO main cycle for current run
		// Keep cycle going until maximum number of evaluations is reached
		do {

			// Update iteration count for current run
			iter += 1;

			// Let particles know about best and worst fitness and determine
			// average fitness
			updatePopulationData(pso);

			// Update all particles
			updateParticles(iter, pso);

			// Is it time to update the average (between runs) best so far?
			if (pso->evaluations > counter * 100) {
				averageBestSoFar[counter] +=
					pso->best_so_far / (long double) n_runs;
				counter += 1;
			}

			// Is the best so far below the stop criteria? If so did we already
			// saved the number of evaluations required to get below the
			// stop criteria?
			if ((pso->best_so_far < crit) && (flag == 0)) {
				crit_evals[i] = pso->evaluations;
				flag = 1;
				// Stop current run if I'm not supposed to keep going
				if (!crit_keep_going) break;
			}

		} while (pso->evaluations < max_evaluations);

		// If the number of evaluations was not enough to get below the stop
		// criteria, set it to the maximum number of performed evaluations
		if (flag == 0) crit_evals[i] = max_evaluations;

		// Inform user of current run performance
		printf("Run %4u | BestFit = %10.5g | AvgFit = %10.5g\n",
			i, (float) pso->best_fitness, (float) pso->average_fitness);

		// Keep best so far for current run
		best_so_far[i] = (double) pso->best_so_far;

		// Release PSO model for current run
		free(pso);
	}

	// Save file containing the average (between runs) best so far fitness
	// after 100, 200, ... evaluations.
	out = fopen("AVE_BESTSOFAR.DAT", "w");
	for (i = 0; i < counter; ++i)
		fprintf(out, "%.40g\n", (double) averageBestSoFar[i]);
	fclose(out);

	// Save number of evaluations required for getting below stop criterion
	out = fopen("AES.DAT", "w");
	for (i = 0; i < n_runs; ++i)
		fprintf(out, "\n%d", crit_evals[i]);
	fclose(out);

	// Save best so far for each run
	out = fopen("FINAL.DAT", "w");
	for (i = 0; i < n_runs; ++i)
		fprintf(out, "%.45g\n", best_so_far[i]);
	fclose(out);

	// Release PRNG states
	free(prng_states);

	// End program successfully
	return EXIT_SUCCESS;
}
