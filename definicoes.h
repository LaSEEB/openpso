/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef __definicoes_h_
#define __definicoes_h_

//#include <assert.h>

#define MAX_POP_SIZE	8100
#define MAX_X			100
#define MAX_Y			100
#define MAX_VARIABLES   100

#define NEIGH_MOORE \
	{{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 0}, \
	{0, 1}, {1, -1}, {1, 0}, {1, 1}}
	//{{0, 0}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, \
	//{0, -1}, {-1, -1}, {-1, 0}, {-1, 1}}
#define NEIGH_MOORE_N 9
#define NEIGH_VN \
	{{0, 0}, {0, 1}, {1, 0}, {0, -1}, {-1, 0}}
#define NEIGH_VN_N 5

typedef struct {
	int dx;
	int dy;
} NEIGHBOR;

typedef struct {
	unsigned int num_neighs;
	const NEIGHBOR * neighs;
} NEIGHBORHOOD;

/////////////////////////////////////////////////////// STRUCTURES
typedef float PARTICLE[MAX_VARIABLES];

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

typedef struct {
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
} MODEL;


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

#endif
