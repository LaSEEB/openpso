/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>

float **aloc_matrizf(int linhas , int colunas);
int **aloc_matrizi(int linhas , int colunas);
int *aloc_vetori(int linhas);
long *aloc_vetorl(int linhas);
float *aloc_vetorf(int linhas);
long double *aloc_vetorld(int linhas);
void desaloc_matrizf(float **Matriz , int linhas);
void desaloc_matrizi(int **Matriz , int linhas);
