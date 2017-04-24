/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "aloca.h"
#include <stdio.h>
#include <stdlib.h>

//These functions allocate memory for data structures.

float ** aloc_matrizf(unsigned int linhas , unsigned int colunas) {
	unsigned int i;
	float ** matrix;
	matrix = (float **) malloc(linhas * sizeof(float));
	for (i = 0; i < linhas; i++) {
		matrix[i] = (float *) malloc(colunas * sizeof(float));
	}
	if (!matrix) {
		printf("Erro de alocacao da matriz ( %d x %d )!", linhas, colunas);
		exit(EXIT_FAILURE);
	}
	return matrix;
}

int ** aloc_matrizi(unsigned int linhas, unsigned int colunas) {
	unsigned int i;
	int ** matrix;
	matrix = (int **)malloc(linhas * sizeof(int));
	for (i = 0; i < linhas; i++) {
		matrix[i] = (int *) malloc(colunas * sizeof(int));
	}
	if (!matrix) {
		printf("Erro de alocacao da matriz ( %d x %d )!", linhas, colunas);
		exit(EXIT_FAILURE);
	}
	return matrix;
}

int * aloc_vetori(unsigned int linhas) {
	int * vetor;
	vetor = (int *) malloc(linhas*sizeof(int));
	if (!vetor) {
		printf("Erro de alocacao de vetor de inteiros!\n");
		exit(EXIT_FAILURE);
	}
	return vetor;
}

long * aloc_vetorl(unsigned int linhas) {
	long * vetor;

	vetor = (long *) malloc(linhas * sizeof(long));
	if (!vetor) {
		printf("Erro de alocacao de vetor de inteiros!\n");
		exit(EXIT_FAILURE);
	}
	return vetor;
}

float *aloc_vetorf(unsigned int linhas) {
	float * vetor;
	vetor = (float *) malloc(linhas * sizeof(float));
	if (!vetor) {
		printf("Erro de alocacao de vetor de floats!\n");
		exit(EXIT_FAILURE);
	}
	return vetor;
}

long double * aloc_vetorld(unsigned int linhas) {
	long double * vetor;
	vetor = (long double *) malloc(linhas * sizeof(long double));
	if (!vetor) {
		printf("Erro de alocacao de vetor de floats!\n");
		exit(EXIT_FAILURE);
	}
	return vetor;
}

void desaloc_matrizf(float ** Matriz , unsigned int linhas) {
	unsigned int i;
	for(i = 0; i < linhas; i++) {
		free(Matriz[i]);
	}
}

void desaloc_matrizi(int ** Matriz , unsigned int linhas) {
	unsigned int i;
	for(i = 0; i < linhas; i++) {
		free(Matriz[i]);
	}
}
