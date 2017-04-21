/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "aloca.h"
#include <stdio.h>
#include <stdlib.h>

//These functions allocate memory for data structures.

float **aloc_matrizf(int linhas , int colunas) {
	int i;
	float **matrix;
	matrix = (float **)malloc(linhas*sizeof(float));
	for (i=0;i<linhas;i++) {
		matrix[i] = (float *)malloc(colunas*sizeof(float));
	}
	if (!matrix) {
		printf("Erro de alocacao da matriz ( %d x %d )!", linhas, colunas);
		exit(1);
	}
	return matrix;
}

int **aloc_matrizi(int linhas , int colunas) {
	int i, **matrix;
	matrix = (int **)malloc(linhas*sizeof(int));
	for (i=0;i<linhas;i++) {
		matrix[i] = (int *)malloc(colunas*sizeof(int));
	}
	if (!matrix) {
		printf("Erro de alocacao da matriz ( %d x %d )!", linhas, colunas);
		exit(1);
	}
	return matrix;
}

int *aloc_vetori(int linhas) {
	int *vetor;
	vetor = (int *)malloc(linhas*sizeof(int));
	if (!vetor) {
		printf("Erro de alocacao de vetor de inteiros!\n");
		exit(1);
	}
	return vetor;
}

long *aloc_vetorl(int linhas) {
	long *vetor;

	vetor = (long *)malloc(linhas*sizeof(long));
	if (!vetor) {
		printf("Erro de alocacao de vetor de inteiros!\n");
		exit(1);
	}
	return vetor;
}

float *aloc_vetorf(int linhas) {
	float *vetor;
	vetor = (float *)malloc(linhas*sizeof(float));
	if (!vetor) {
		printf("Erro de alocacao de vetor de floats!\n");
		exit(1);
	}
	return vetor;
}

long double *aloc_vetorld(int linhas) {
	long double *vetor;
	vetor = (long double *)malloc(linhas*sizeof(long double));
	if (!vetor) {
		printf("Erro de alocacao de vetor de floats!\n");
		exit(1);
	}
	return vetor;
}

void desaloc_matrizf(float **Matriz , int linhas) {
	int i;
	for(i=0;i<linhas;i++) {
		free(Matriz[i]);
	}
}

void desaloc_matrizi(int **Matriz , int linhas) {
	int i;
	for(i = 0; i < linhas; i++) {
		free(Matriz[i]);
	}
}
