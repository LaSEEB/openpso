/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

///////////////////////////////////////////////////////
// SS-PSO
// Carlos Fernandes
// Nuno Fachada

//#include <malloc.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include "definicoes.h"
#include "aloca.h"
#include "functions.h"

// *************** Global Variables *****************

// System's parameters

int max_x, max_y;
// habitat x*y
int max_t;
// Maximum number of iterations
unsigned int max_evaluations;
// Population size
int popSize;
// Number of particles in the swarm at time t=0
unsigned int n_runs;
// Array used to check the occupation of each ant' surroundig sites
// ? Run Bak-Sneppen model on the structure
int position[8][2];
// run PSO on the structure
int algorithm;
//
int gbest;
//
int neighborhood;
//
int problem;
//
float Xmax;
//
float Vmax;
//
float chi;
//
float omega;
//
float c;
//
int numberVariables;
//
int assyInitialization;
//
float initialXmin;
//
float initialXmax;
//
float crit;
//
int iWeightStrategy,cStrategy;
//
int parametros, semilla, pausa;


/////RANDOM NUMBERS AND MATH FUNCTIONS
////////////////////////////////////////////////////////////
long double	real_al_entre_a_b(long double a,long double b) {
	return  a+(rand()*(b-a)) / (pow(2.0,15.0)-1.0);
}
////////////////////////////////////////////////////////////
long aleatorio_entre_a_b(long double a,long double b) {
	return	(long) floor(a+(rand()*(b+1-a)) / (pow(2.0,15.0)));
}
////////////////////////////////////////////////////////////
double square(float x) {
	return (x*x);
}
/////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// PARTICLE SWARM OPTIMIZATION    ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
void initialize(MODEL * pso) {

	int i,j, z;
	float xmin, xmax;

	z = 0;
	for (i = 0; i < max_x; ++i){
		for (j = 0; j < max_y; ++j){
			pso->cell[i][j].particle = z;
			pso->particle[z].x = i;
			pso->particle[z].y = j;
			z = z+1;
		}
	}

	//Initialize Position and Velocity
	for (i = 0;  i < popSize;    ++i) {

		if (assyInitialization == 1) {
			// Assymetric initialization of the population
		   xmin = initialXmin;
		   xmax = initialXmax;
		} else {
			// Normal initialization
			 xmin = -Xmax;
			 xmax = Xmax;
		}

		for (j = 0;  j < numberVariables;   ++j) {
			pso->particle[i].position[j] = real_al_entre_a_b(xmin, xmax);
			pso->particle[i].velocity[j] = real_al_entre_a_b(-Xmax, Xmax)*(0.5-real_al_entre_a_b(0,1.0));
			pso->particle[i].best_position_so_far[j] = pso->particle[i].position[j];
			pso->particle[i].informants_best_position_so_far[j] = pso->particle[i].position[j];
		}

		pso->particle[i].fitness = evaluate(pso->particle[i].position);
		pso->particle[i].best_fitness_so_far = pso->particle[i].fitness;
		pso->particle[i].informants_best_fitness_so_far = pso->particle[i].fitness;
	}
	pso->best_so_far = pso->particle[0].fitness;
	pso->best_so_far_id = 0;
	pso->evaluations = 0;
	pso->minFitness = pso->particle[0].fitness;

}

///////////////////////////////////
void updatePopulationData (MODEL *pso) {

	int i,j;
	pso->best_fitness = pso->particle[0].fitness;
	pso->worst_fitness = 0;
	pso->average_fitness = 0;
	for (i = 0;  i < popSize;     ++i) {
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
			   pso->particle[i].best_position_so_far[j] = pso->particle[i].position[j];
		}
		// Updates best informant
		if (pso->particle[i].fitness < pso->particle[i].informants_best_fitness_so_far) {
		   pso->particle[i].informants_best_fitness_so_far = pso->particle[i].fitness;
		   for (j = 0;  j < numberVariables;   ++j)
			   pso->particle[i].informants_best_position_so_far[j] = pso->particle[i].position[j];
		}
		pso->average_fitness = pso->average_fitness+pso->particle[i].fitness;
	}
	pso->average_fitness = pso->average_fitness/popSize;
}

/////////////////////////////////////
void updateParticles (MODEL *pso, int i, unsigned int t) {
	int j;
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
		phi1 = real_al_entre_a_b (0.0, c1);
		phi2 = real_al_entre_a_b (0.0, c2);
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

	int i, a, j, ii, jj, z;
	int minx, maxx, miny, maxy;
	int neighborAnt;
	int update;
	 for (a = 0;   a < popSize;       ++a) {
		 update = 0;
		 minx = pso->particle[a].x-1;
		 miny = pso->particle[a].y-1;
		 maxx = pso->particle[a].x+1;
		 maxy = pso->particle[a].y+1;
		 for (i = minx;	i <= maxx;	++i) {
			 for (j = miny;	j <= maxy;	++j) {
				 ii = i;
				 jj = j;
				 if (i < 0)      ii = max_x-1;
				 if (i >= max_x) ii = 0;
				 if (j < 0)      jj = max_y-1;
				 if (j >= max_y) jj = 0;
				 // Updates best neighbor
				 if (neighborhood == 0 || (i == minx+1 && j == miny+1) || (i == minx+1 && j == miny) || (i == minx && j == miny+1) || (i == minx+1 && j == maxy) || (i == maxx && j == miny+1)) {
					if (pso->cell[ii][jj].particle == pso->worst_id) // the worst ant is a neighbor
						update = 1;    // mark particle for updating; used for SS-PSO
					neighborAnt = pso->cell[ii][jj].particle;
					if (pso->particle[neighborAnt].best_fitness_so_far < pso->particle[a].informants_best_fitness_so_far) {
						pso->particle[a].informants_best_fitness_so_far = pso->particle[neighborAnt].best_fitness_so_far;
						for (z = 0;   z < numberVariables; ++z)
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

//################################################################## MAIN ALGORITHM
int main(int argc, char* argv[]) {

	FILE *out1;
	unsigned int i, w;
	MODEL *pso;
	unsigned int z;
	unsigned int counter = 0;
	FILE *in1;
	long double *averageBestsofar;
	int flag;
	// Read parameters ////////////
	in1=fopen("INPUT.txt","r");
	fscanf(in1,"%d", &max_x);
	fscanf(in1,"%d", &max_y);
	fscanf(in1,"%d", &n_runs);
	fscanf(in1,"%d", &max_t);
	fscanf(in1,"%d", &max_evaluations);
	fscanf(in1,"%d", &popSize);
	fscanf(in1,"%d", &algorithm);
	fscanf(in1,"%d", &gbest);
	fscanf(in1,"%d", &neighborhood);
	fscanf(in1,"%d", &problem);
	fscanf(in1,"%f", &Xmax);
	fscanf(in1,"%f", &Vmax);
	fscanf(in1,"%f", &chi);
	fscanf(in1,"%f", &omega);
	fscanf(in1,"%f", &c);
	fscanf(in1,"%u", &numberVariables);
	fscanf(in1,"%d", &iWeightStrategy);
	fscanf(in1,"%d", &cStrategy);
	fscanf(in1,"%d", &assyInitialization);
	fscanf(in1,"%f", &initialXmin);
	fscanf(in1,"%f", &initialXmax);
	fscanf(in1,"%f", &crit);
	fclose (in1);
	//////////////////////////////
	averageBestsofar = aloc_vetorld(50000);
	for (i = 0;      i < 50000;           ++i)
		averageBestsofar[i] = 0.0;
	for (w = 0;   w < n_runs;   ++w) {
	   if((pso = (MODEL *) calloc(1,sizeof(MODEL))) == NULL) {
			printf("\nERROR: Out of Memory");
			exit(0);
		}
		flag = 0;
		z = 0;
		counter = 0;
		initialize(pso);
		// *************** MAIN CYCLE *****************
		do {
			z = z+1;
			updatePopulationData (pso);
			printf ("\n best_so_far = %f", (float)pso->average_fitness);
			move(z, pso);
			if (pso->evaluations > counter * 100) {
			   averageBestsofar[counter] = averageBestsofar[counter]+pso->best_so_far/(long double)n_runs;
			   counter = counter+1;
			}
			if (pso->best_so_far < crit && flag == 0) {
				out1=fopen("AES.DAT","a");
				fprintf(out1,"\n%d", pso->evaluations);
				fclose (out1);
				flag = 1;
			}
		} while (pso->evaluations < max_evaluations);
		out1=fopen("INTERMEDIARY.DAT","a");
		fprintf(out1,"\n");
		fclose (out1);
		out1=fopen("FINAL.DAT","a");
		fprintf(out1,"%.45f\n", (float)pso->best_so_far);
		fclose (out1);
		//getch();
		free (pso);
	}
	out1=fopen("AVE_BESTSOFAR.DAT","a");
	for (i = 1;	i < counter+1;	++i)
	   fprintf(out1,"%.40f\n", (float)averageBestsofar[i]);
	fclose (out1);
	free(averageBestsofar);
	return 0;
} // main
/////////////////////////////////////////////////////// END PSO ALGORITHM
