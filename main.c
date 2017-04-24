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
#include "aloca.h"
#include "functions.h"

// Constants
#define DEFAULT_INPUT_FILE "input.ini"
#define DEFAULT_PRNG_SEED 1234

// *************** Global Variables *****************

/// Seed for pseudo-random number generator.
static unsigned int prng_seed;
/// File containing PSO parameters.
static char * input_file;

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
unsigned int problem;
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
unsigned int numberVariables;
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
 * Macro for terminating program with error condition.
 *
 * @param[in] err_msg Error message.
 * @param[in] ... Message parameters.
 */
#define ERROR_EXIT(err_msg, ...) \
	do { \
		fprintf(stderr, err_msg, ##__VA_ARGS__); \
		exit(EXIT_FAILURE); \
	} while(0)


//////////////////////////////////////////////////////////////////////////////////
/////////////////////////// PARTICLE SWARM OPTIMIZATION    ///////////////////////
//////////////////////////////////////////////////////////////////////////////////
void initialize(MODEL * pso) {

	unsigned int i, j, z;
	float xmin, xmax;

	z = 0;
	for (i = 0; i < max_x; ++i){
		for (j = 0; j < max_y; ++j){
			pso->cell[i][j].particle = z;
			pso->particle[z].x = i;
			pso->particle[z].y = j;
			z = z + 1;
		}
	}

	//Initialize Position and Velocity
	for (i = 0; i < popSize; ++i) {

		if (assyInitialization == 1) {
			// Assymetric initialization of the population
			xmin = initialXmin;
			xmax = initialXmax;
		} else {
			// Normal initialization
			xmin = -Xmax;
			xmax = Xmax;
		}

		for (j = 0;  j < numberVariables; ++j) {
			pso->particle[i].position[j] =
				rd_luniform(xmin, xmax);
			pso->particle[i].velocity[j] =
				rd_luniform(-Xmax, Xmax)*(0.5-rd_luniform(0,1.0));
			pso->particle[i].best_position_so_far[j] =
				pso->particle[i].position[j];
			pso->particle[i].informants_best_position_so_far[j] =
				pso->particle[i].position[j];
		}

		pso->particle[i].fitness = evaluate(pso->particle[i].position);
		pso->particle[i].best_fitness_so_far = pso->particle[i].fitness;
		pso->particle[i].informants_best_fitness_so_far = pso->particle[i].fitness;
	}
	pso->best_so_far = pso->particle[0].fitness;
	pso->best_so_far_id = 0;
	pso->evaluations = 0;
	pso->minFitness = pso->particle[0].fitness;
	pso->worst_id = 0;

}

///////////////////////////////////
void updatePopulationData (MODEL *pso) {

	unsigned int i,j;
	pso->best_fitness = pso->particle[0].fitness;
	pso->worst_fitness = 0;
	pso->average_fitness = 0;
	for (i = 0; i < popSize; ++i) {
		// Updates worst in population
		if (pso->particle[i].fitness > pso->worst_fitness) {
			pso->worst_fitness = pso->particle[i].fitness;
			pso->worst_id = i;
		}
		// Updates best_so_far in population
		if (pso->particle[i].fitness < pso->best_so_far) {
			pso->best_so_far = pso->particle[i].fitness;
			for (j = 0;   j < numberVariables;   ++j)
				pso->best_position_so_far[j] = pso->particle[i].position[j];
		}
		// Updates best in current population
		if (pso->particle[i].fitness < pso->best_fitness) {
			pso->best_fitness = pso->particle[i].fitness;
			for (j = 0;   j < numberVariables;   ++j)
				pso->best_position[j] = pso->particle[i].position[j];
		}
		// Updates particle's best position
		if (pso->particle[i].fitness < pso->particle[i].best_fitness_so_far) {
			pso->particle[i].best_fitness_so_far = pso->particle[i].fitness;
			for (j = 0;  j < numberVariables;  ++j)
				pso->particle[i].best_position_so_far[j] =
				pso->particle[i].position[j];
		}
		// Updates best informant
		if (pso->particle[i].fitness <
				pso->particle[i].informants_best_fitness_so_far) {

			pso->particle[i].informants_best_fitness_so_far =
				pso->particle[i].fitness;
			for (j = 0; j < numberVariables; ++j)
				pso->particle[i].informants_best_position_so_far[j] =
					pso->particle[i].position[j];

		}
		pso->average_fitness = pso->average_fitness + pso->particle[i].fitness;
	}
	pso->average_fitness = pso->average_fitness / popSize;
}

/////////////////////////////////////
void updateParticles (MODEL *pso, int i, unsigned int t) {
	unsigned int j;
	float v, x;
	float pi, pg;
	long double phi1, phi2;
	float c1, c2;
	float maxIW, minIW;
	FILE *out1;
	c1 = c;
	c2 = c;
	maxIW = 0.9;
	minIW = 0.4;
	if (iWeightStrategy == 1)     // TVIW-PSO
		omega = minIW+(maxIW-minIW)*(((float)max_t-(float)t)/(float)max_t);
	if (cStrategy == 1) {         //TVAC-PSO
		c1 = (0.5-2.5)*((float)t/(float)max_t)+2.5;
		c2 = (2.5-0.5)*((float)t/(float)max_t)+0.5;
	}
	for (j = 0;  j < numberVariables;    ++j) {
		if (gbest == 0) pg = pso->particle[i].informants_best_position_so_far[j];
		if (gbest == 1) pg = pso->best_position_so_far[j];
		pi = pso->particle[i].best_position_so_far[j];
		v = pso->particle[i].velocity[j];
		x = pso->particle[i].position[j];
		phi1 = rd_luniform(0.0, c1);
		phi2 = rd_luniform(0.0, c2);
		// Update Velocity
		v = omega*v+(float)phi1*(pi-x)+(float)phi2*(pg-x);
		if (v > Vmax) v = Vmax;
		if (v < -Vmax) v = -Vmax;
		// Update Position
		x = x+v;
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
	pso->particle[i].fitness = evaluate (pso->particle[i].position);
	pso->evaluations = pso->evaluations+1;
	if (pso->evaluations == 49000 || pso->evaluations == 147000 || pso->evaluations == 294000 || pso->evaluations == 490000) {
		out1=fopen("INTERMEDIARY.DAT","a");
		fprintf(out1,"%.50f\t", (float)pso->best_so_far);
		fclose (out1);
	}
}
// END PARTICLE SWARM OPTIMIZATION

//////////////////////MOVE THE PARTICLES
void move(unsigned int t, MODEL *pso) {

	int a, i, j, ii, jj;
	unsigned int z;
	int minx, maxx, miny, maxy;
	int neighborAnt;
	int update;
	for (a = 0; a < (int) popSize; ++a) {
		update = 0;
		minx = pso->particle[a].x-1;
		miny = pso->particle[a].y-1;
		maxx = pso->particle[a].x+1;
		maxy = pso->particle[a].y+1;
		for (i = minx; i <= maxx; ++i) {
			for (j = miny; j <= maxy; ++j) {
				ii = i;
				jj = j;
				if (i < 0)      ii = max_x - 1;
				if (i >= (int) max_x) ii = 0;
				if (j < 0)      jj = max_y - 1;
				if (j >= (int) max_y) jj = 0;
				// Updates best neighbor
				if (neighborhood == 0 || (i == minx+1 && j == miny+1) || (i == minx+1 && j == miny) || (i == minx && j == miny+1) || (i == minx+1 && j == maxy) || (i == maxx && j == miny+1)) {
					if (pso->cell[ii][jj].particle == pso->worst_id) // the worst ant is a neighbor
						update = 1;    // mark particle for updating; used for SS-PSO
					neighborAnt = pso->cell[ii][jj].particle;
					if (pso->particle[neighborAnt].best_fitness_so_far < pso->particle[a].informants_best_fitness_so_far) {
						pso->particle[a].informants_best_fitness_so_far = pso->particle[neighborAnt].best_fitness_so_far;
						for (z = 0; z < numberVariables; ++z)
							pso->particle[a].informants_best_position_so_far[z] = pso->particle[neighborAnt].best_position_so_far[z];
					}
				}
		   }
		}
		if (algorithm == 1)
			updateParticles (pso, a, t);
		if (algorithm == 2 && update == 1) // new PSO; only the worst and its neigbobors are updated
			updateParticles (pso, a, t);
	 }
}

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
		input_file = argv[1];
	else
		input_file = DEFAULT_INPUT_FILE;

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
 * Algorithm starts here.
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

	averageBestSoFar = aloc_vetorld(50000);
	for (i = 0; i < 50000; ++i)
		averageBestSoFar[i] = 0.0;

	for (i = 0; i < n_runs; ++i) {

		// Create PSO model, set contents to zero
		pso = (MODEL *) calloc(1, sizeof(MODEL));

		flag = 0;
		z = 0;
		counter = 0;
		initialize(pso);

		// *************** MAIN CYCLE *****************
		do {
			z = z + 1;
			updatePopulationData(pso);
			printf ("\n best_so_far = %f", (float) pso->average_fitness);
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

		//getch();
		free(pso);
	}

	out1 = fopen("AVE_BESTSOFAR.DAT", "a");
	for (i = 1; i < counter + 1; ++i)
		fprintf(out1,"%.40f\n", (float) averageBestSoFar[i]);
	fclose (out1);

	free(averageBestSoFar);

	return 0;
}
