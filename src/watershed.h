/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Definitions of for the watershed update strategy.
 *
 * @author Nuno Fachada
 */

#ifndef __PSO_WATERSHED_H_
#define __PSO_WATERSHED_H_

#include "common.h"

/// Watershed minimum fitness strategy function
typedef void (*pso_watershed)(PSO *pso);

/// Watershed strategy selection
void pso_watershed_select(PSO *pso, const char *strategy);

/// Existing functions
void pso_watershed_worst_so_far(PSO *pso);
void pso_watershed_worst_last_iter(PSO *pso);
void pso_watershed_best_worst(PSO *pso);

#endif
