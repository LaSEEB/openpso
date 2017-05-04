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
#include <limits.h>
#ifdef _OPENMP
	#include <omp.h>
#endif

// Local libraries
#include "randistrs.h"
#include "mtwist.h"
#include "iniparser.h"

// SSPSO headers
#include "functions.h"
#include "pso.h"

// Constants
#define DEFAULT_INPUT_FILE "input.ini"
#define DEFAULT_PRNG_SEED 1234
#define MAX_FILENAME_LEN 255
/// Maximum length for strings containing unsigned integers
#define MAXUISLEN 11

// Useful macros

/**
 * Macro for terminating program with error condition.
 *
 * @param[in] err_msg Error message.
 * @param[in] ... Message parameters.
 */
#define ERROR_EXIT(err_msg, ...) \
	do { \
		fprintf(stderr, err_msg, ##__VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while (0)

/**
 * Macro for calculating median. Assumes sorted vector.
 *
 * @param[in] vec Sorted numeric vector.
 * @param[in] n Number of elements in vector.
 * @param[out] res Median of numeric elements in vector.
 */
#define MEDIAN(vec, n, res) \
	do { \
		if (n % 2 == 0) res = (vec[n / 2] + vec[n / 2 - 1]) / 2; \
		else res = vec[n / 2]; \
	} while (0)

static pso_func getSelFunc(unsigned int func) {
	switch (func) {
		case 1: return &Sphere;
		case 2: return &Quadric;
		case 3: return &Hyper;
		case 4: return &Rastrigin;
		case 5: return &Griewank;
		case 6: return &Schaffer6;
		case 7: return &Weierstrass;
		case 8: return &Ackley;
		case 9: return &ShiftedQuadricWithNoise;
		case 10: return &RotatedGriewank;
		default: return NULL;
	}
}

// *************** Global Variables *****************



/// Seed for pseudo-random number generator.
static unsigned int prng_seed;
/// File containing PSO parameters.
static char input_file[MAX_FILENAME_LEN];

/// Number of threads.
static unsigned int num_threads;
/// Problem to solve
unsigned int problem;
// Number of runs
unsigned int n_runs;
// Maximum number of evaluations
unsigned int max_evaluations;

// Helper function for comparing two doubles
static int cmpdbl(const void * a, const void * b) {
	if (*(double*) a > *(double*) b) return 1;
	else if (*(double*) a < *(double*) b) return -1;
	else return 0;
}

// Helper function for comparing two unsigned integers
static int cmpuint(const void * a, const void * b) {
	if (*(unsigned int*) a > *(unsigned int*) b) return 1;
	else if (*(unsigned int*) a < *(unsigned int*) b) return -1;
	else return 0;
}

// Helper function for converting unsigned integer to string
// Cannot be called by threaded code
static char * uint2str(unsigned int evals) {
	static char str[MAXUISLEN]; // Not thread-safe
	snprintf(str, MAXUISLEN * sizeof(char), "%u", evals);
	return str;
}

/**
 * Parse command-line options and read PSO parameters from file.
 *
 * @param[in] argc Number of program arguments.
 * @param[in] argv Program arguments. Argument at index 1 is the file containing
 * the PSO parameters. Argument at index 2 is the PRNG seed.
 * @param[out] params Pointer to PSO parameters.
 * @return `EXIT_SUCCESS` if program executes successfully, `EXIT_FAILURE`
 * otherwise.
 */
static void parse_params(int argc, char * argv[], PARAMETERS * params) {

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

	// Try to open PSO parameters file
	ini = iniparser_load(input_file);
	if (!ini) ERROR_EXIT("Unable to parse input file '%s'", input_file);

	// Read PSO parameters file

	params->max_x = (unsigned int) iniparser_getint(ini, "pso:max_x", 0);

	params->max_y = (unsigned int) iniparser_getint(ini, "pso:max_y", 0);

	n_runs = (unsigned int) iniparser_getint(ini, "pso:n_runs", 0);
	if (n_runs < 1)
		ERROR_EXIT("Invalid input parameter: n_runs");

	params->max_t = (unsigned int) iniparser_getint(ini, "pso:max_t", 0);

	max_evaluations = (unsigned int) iniparser_getint(ini, "pso:max_evaluations", 0);
	if (max_evaluations < 1)
		ERROR_EXIT("Invalid input parameter: max_evaluations");

	params->algorithm = (unsigned int) iniparser_getint(ini, "pso:algorithm", 0);

	params->gbest = iniparser_getboolean(ini, "pso:gbest", -1);

	// Neighborhood
	params->neighborhood = (unsigned int) iniparser_getint(ini, "pso:neighborhood", 3);

	problem = (unsigned int) iniparser_getint(ini, "pso:problem", 0);
	if ((problem < 1) || (problem > 10))
		ERROR_EXIT("Invalid input parameter: problem");

	params->Xmax = iniparser_getdouble(ini, "pso:xmax", -DBL_MAX);

	params->Vmax = iniparser_getdouble(ini, "pso:vmax", -DBL_MAX);

	params->chi = iniparser_getdouble(ini, "pso:chi", -DBL_MAX);

	params->omega = iniparser_getdouble(ini, "pso:omega", -DBL_MAX);

	params->c = iniparser_getdouble(ini, "pso:c", -DBL_MAX);

	params->nvars = (unsigned int) iniparser_getint(ini, "pso:numbervariables", 0);

	params->iWeightStrategy = (unsigned int) iniparser_getint(ini, "pso:iweightstrategy", 2);

	params->cStrategy = (unsigned int) iniparser_getint(ini, "pso:cstrategy", 2);

	params->assyInitialization = iniparser_getboolean(ini, "pso:assyinitialization", -1);

	params->initialXmin = iniparser_getdouble(ini, "pso:initialxmin", -DBL_MAX);

	params->initialXmax = iniparser_getdouble(ini, "pso:initialxmax", -DBL_MAX);

	params->crit = iniparser_getdouble(ini, "pso:crit", -DBL_MAX);

 	params->crit_keep_going = iniparser_getboolean(ini, "pso:crit_keep_going", -1);

	// Release dictionary object
	iniparser_freedict(ini);

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

	PARAMETERS params;

	// Parse command-line arguments and read PSO parameter file.
	parse_params(argc, argv, &params);

	FILE * out;
	unsigned int i;
	PSO * pso;

	// Iteration count for current run
	unsigned int iter;

	// Medians for fitness and evaluations
	float fitmedian, evalsmedian;

	unsigned int counter = 0;
	long double averageBestSoFar[max_evaluations / 100];
	unsigned int crit_evals[n_runs];
	unsigned int successes = 0;
	double best_so_far[n_runs];
	int flag;

	// Show parameter information to user
	printf("\nPARAMS\n------\n");
	printf("Num threads        : %u\n", num_threads);
	printf("PSO parameter file : %s\n", input_file);
	printf("PRNG seed          : %d\n", prng_seed);
	printf("\n");

	// Initialize averageBestSoFar array, set contents to zero
	memset(averageBestSoFar, 0, max_evaluations / 100 * sizeof(long double));

	printf("\nRUNS\n----\n");

	// Perform PSO runs
	for (i = 0; i < n_runs; ++i) {

		// Initialize PSO for current run
		pso = pso_new(params , getSelFunc(problem), prng_seed);

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
			if ((pso->best_so_far < pso->params.crit) && (flag == 0)) {
				crit_evals[i] = pso->evaluations;
				flag = 1;
				successes++;
				// Stop current run if I'm not supposed to keep going
				if (!pso->params.crit_keep_going) break;
			}

		} while (pso->evaluations < max_evaluations);

		// If the number of evaluations was not enough to get below the stop
		// criteria, set it to the maximum number of performed evaluations
		if (flag == 0) {
			crit_evals[i] = UINT_MAX;
		}

		// Inform user of current run performance
		printf("Run %4u | BestFit = %10.5g | AvgFit = %10.5g | "
			"Evals = %10s\n", i, (double) pso->best_so_far,
			(double) pso->average_fitness,
			flag ? uint2str(crit_evals[i]) : "--");

		// Keep best so far for current run
		best_so_far[i] = (double) pso->best_so_far;

		// Release PSO model for current run
		pso_destroy(pso);
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
	for (i = 0; i < n_runs; ++i) {
		fprintf(out, "%.45g\n", best_so_far[i]);
	}
	fclose(out);

	// Show statistics over all runs
	printf("\nSTATISTICS\n----------\n");
	printf("          \t %11s \t %11s \t %11s\n", "Median", "Min", "Max");

	// Fitness statistics
	qsort(best_so_far, n_runs, sizeof(double), cmpdbl);
	MEDIAN(best_so_far, n_runs, fitmedian);

	// Evaluation statistics
	printf("Fitness = \t %10.5e \t %10.5e \t %10.5e\n",
		fitmedian, best_so_far[0], best_so_far[n_runs - 1]);
	if (successes) {
		qsort(crit_evals, n_runs, sizeof(unsigned int), cmpuint);
		MEDIAN(crit_evals, successes, evalsmedian);
		printf("Evals   = \t %11.1f \t %11u \t %11u\n",
			evalsmedian, crit_evals[0], crit_evals[successes - 1]);
	} else {
		printf("Evals   = \t %11s \t %11s \t %11s\n", "--", "--", "--");
	}

	// Number and percentage of successes
	printf("\nSuccess = %u/%u (%2.2f%%)\n\n", successes, n_runs,
		100.0 * (float) successes / (float) n_runs);

	// End program successfully
	return EXIT_SUCCESS;
}
