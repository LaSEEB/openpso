/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "definicoes.h"

void initialize(MODEL *kants);
void leer_patrones(char *nombre_fich, MODEL *kants);
void normalizar_valores(MODEL *kants);
