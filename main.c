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

/// Seed for pseudo-random number generator.
static unsigned int prng_seed;
/// File containing PSO parameters.
static char input_file[MAX_FILENAME_LEN];
// Function/problem to solve
static SelFunc evaluate;

// PSO parameters

//
static unsigned int max_x, max_y;
//Maximum number of iterations
static unsigned int max_t;
// Maximum number of evaluations
static unsigned int max_evaluations;
// Population size
static unsigned int popSize;
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
//
static double crit;

/**
 * Parse command-line options and read PSO parameters from file.
 *
 * @param[in] argc Number of program arguments.
 * @param[in] argv Program arguments. Argument at index 1 is the file containing
 * the PSO parameters. Argument at index 2 is the PRNG seed.
 * @return `EXIT_SUCCESS` if program executes successfully, `EXIT_FAILURE`
 * otherwise.
 */
void parse_params(int argc, char * argv[]) {

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
	popSize = (unsigned int) iniparser_getint(ini, "pso:popsize", 0);
	if (popSize < 1)
		ERROR_EXIT("Invalid input parameter: popSize\n");
	algorithm = (unsigned int) iniparser_getint(ini, "pso:algorithm", 0);
	if ((algorithm < 1) || (algorithm > 2))
		ERROR_EXIT("Invalid input parameter: algorithm\n");
	gbest = iniparser_getboolean(ini, "pso:gbest", -1);
	if (gbest == -1)
		ERROR_EXIT("Invalid input parameter: gbest\n");
	neighborhood = (unsigned int) iniparser_getint(ini, "pso:neighborhood", 2);
	if (neighborhood == 2)
		ERROR_EXIT("Invalid input parameter: neighborhood\n");
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

	// Release dictionary object
	iniparser_freedict(ini);

}

/**
 * Initialize PSO model.
 *
 * @param[in,out] PSO model to initialize.
 */
void initialize(MODEL * pso) {

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
				rd_luniform(xmin, xmax);

			// Initialize velocity for current variable of current particle
			pso->particle[i].velocity[j] =
				rd_luniform(-Xmax, Xmax) * (0.5 - rd_luniform(0, 1.0));

			// Set best position so far as current position
			pso->particle[i].best_position_so_far[j] =
				pso->particle[i].position[j];

			// Set best informat so far as myself
			pso->particle[i].informants_best_position_so_far[j] =
				pso->particle[i].position[j];
		}

		// Determine fitness for current particle
		pso->particle[i].fitness =
			evaluate(pso->particle[i].position, numberVariables);

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
void updatePopulationData(MODEL * pso) {

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

/////////////////////////////////////
void updateParticles(MODEL * pso, int i, unsigned int t) {
	unsigned int j;
	float v, x;
	float pi, pg = 0.0;
	long double phi1, phi2;
	float c1, c2;
	float maxIW, minIW;
	FILE * out1;
	c1 = c;
	c2 = c;
	maxIW = 0.9;
	minIW = 0.4;
	if (iWeightStrategy == 1)     // TVIW-PSO
		omega =
			minIW + (maxIW - minIW) * (((float)max_t-(float)t)/(float)max_t);
	if (cStrategy == 1) {         //TVAC-PSO
		c1 = (0.5 - 2.5) * ((float) t / (float) max_t) + 2.5;
		c2 = (2.5 - 0.5) * ((float) t / (float) max_t) + 0.5;
	}
	for (j = 0; j < numberVariables; ++j) {

		if (gbest == 0)
			pg = pso->particle[i].informants_best_position_so_far[j];
		if (gbest == 1)
			pg = pso->best_position_so_far[j];

		pi = pso->particle[i].best_position_so_far[j];
		v = pso->particle[i].velocity[j];
		x = pso->particle[i].position[j];
		phi1 = rd_luniform(0.0, c1);
		phi2 = rd_luniform(0.0, c2);

		// Update Velocity
		v = omega * v + (float) phi1 * (pi - x) + (float) phi2 * (pg - x);
		if (v > Vmax) v = Vmax;
		if (v < -Vmax) v = -Vmax;

		// Update Position
		x = x + v;
		if (x > Xmax) {
			x = Xmax;
			v = 0;
		}
		if (x < -Xmax) {
			x = -Xmax;
			v = 0;
		}
		pso->particle[i].position[j] = x;
		pso->particle[i].velocity[j] = v;
	}
	pso->particle[i].fitness =
		evaluate(pso->particle[i].position, numberVariables);
	pso->evaluations = pso->evaluations + 1;
	if (pso->evaluations == 49000 ||
			pso->evaluations == 147000 ||
			pso->evaluations == 294000 ||
			pso->evaluations == 490000) {

		out1 = fopen("INTERMEDIARY.DAT", "a");
		fprintf(out1, "%.50f\t", (float) pso->best_so_far);
		fclose (out1);

	}
}
// END PARTICLE SWARM OPTIMIZATION

//////////////////////MOVE THE PARTICLES
void move(unsigned int t, MODEL *pso) {

	int a, i, j, ii, jj;
	unsigned int z;
	int minx, maxx, miny, maxy;
	int neighParticle;
	int update;
	for (a = 0; a < (int) popSize; ++a) {
		update = 0;
		minx = pso->particle[a].x - 1;
		miny = pso->particle[a].y - 1;
		maxx = pso->particle[a].x + 1;
		maxy = pso->particle[a].y + 1;

		for (i = minx; i <= maxx; ++i) {

			for (j = miny; j <= maxy; ++j) {

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

				// Updates best neighbor
				if (neighborhood == 0 ||
						(i == minx + 1 && j == miny + 1) ||
						(i == minx + 1 && j == miny) ||
						(i == minx && j == miny + 1) ||
						(i == minx + 1 && j == maxy) ||
						(i == maxx && j == miny + 1)) {

					if (pso->cell[ii][jj].particle == pso->worst_id)
						// the worst particle is a neighbor
						update = 1; // mark particle for updating (for SS-PSO)

					neighParticle = pso->cell[ii][jj].particle;
					if (pso->particle[neighParticle].best_fitness_so_far <
							pso->particle[a].informants_best_fitness_so_far) {
						pso->particle[a].informants_best_fitness_so_far =
							pso->particle[neighParticle].best_fitness_so_far;
						for (z = 0; z < numberVariables; ++z)
							pso->particle[a].informants_best_position_so_far[z] =
								pso->particle[neighParticle].best_position_so_far[z];
					}
				}
		   }
		}

		if (algorithm == 1)
			updateParticles(pso, a, t);

		if ((algorithm == 2) && (update == 1))
			// new PSO; only the worst and its neighbors are updated
			updateParticles (pso, a, t);
	 }
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

	FILE * out1;
	unsigned int i;
	MODEL * pso;
	unsigned int z;
	unsigned int counter = 0;
	long double * averageBestSoFar;
	int flag;

	// Parse command-line arguments and read PSO parameter file.
	parse_params(argc, argv);

	// Show parameter information to user
	printf("PSO parameter file : %s\n", input_file);
	printf("PRNG seed          : %d\n", prng_seed);

	// Set PRNG seed.
	mt_seed32new(prng_seed);

	// Set function/problem to optimize
	evaluate = getSelFunc(problem);

	// Initialize averageBestSoFar array, set contents to zero
	averageBestSoFar = (long double *) calloc(50000, sizeof(long double));

	// Perform PSO runs
	for (i = 0; i < n_runs; ++i) {

		// Create PSO model for current run, set contents to zero
		pso = (MODEL *) calloc(1, sizeof(MODEL));

		// Initialize PSO for current run
		initialize(pso);

		// Aux. variables for current run
		flag = 0;
		z = 0;
		counter = 0;

		// PSO main cycle for current run
		// Keep cycle going until maximum number of evaluations is reached
		do {

			z = z + 1;

			updatePopulationData(pso);

			printf("\n best_so_far = %f", (float) pso->average_fitness);

			move(z, pso);

			if (pso->evaluations > counter * 100) {
				averageBestSoFar[counter] =
					averageBestSoFar[counter] +
					pso->best_so_far / (long double) n_runs;
				counter += 1;
			}

			if ((pso->best_so_far < crit) && (flag == 0)) {
				out1 = fopen("AES.DAT", "a");
				fprintf(out1, "\n%d", pso->evaluations);
				fclose(out1);
				flag = 1;
			}

		} while (pso->evaluations < max_evaluations);

		out1 = fopen("INTERMEDIARY.DAT", "a");
		fprintf(out1, "\n");
		fclose (out1);

		out1 = fopen("FINAL.DAT", "a");
		fprintf(out1, "%.45f\n", (float) pso->best_so_far);
		fclose(out1);

		// Release PSO model for current run
		free(pso);
	}

	out1 = fopen("AVE_BESTSOFAR.DAT", "a");
	for (i = 1; i < counter + 1; ++i)
		fprintf(out1,"%.40f\n", (float) averageBestSoFar[i]);
	fclose (out1);

	free(averageBestSoFar);

	return 0;
}
