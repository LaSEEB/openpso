/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>

float ** aloc_matrizf(unsigned int linhas , unsigned int colunas);
int ** aloc_matrizi(unsigned int linhas , unsigned int colunas);
int * aloc_vetori(unsigned int linhas);
long * aloc_vetorl(unsigned int linhas);
float * aloc_vetorf(unsigned int linhas);
long double * aloc_vetorld(unsigned int linhas);
void desaloc_matrizf(float ** Matriz , unsigned int linhas);
void desaloc_matrizi(int ** Matriz , unsigned int linhas);
