/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Watershed strategy functions.
 *
 * @author Nuno Fachada
 */

#include "watershed.h"
#include "pso.h"
#include "zf_log.h"
#include <string.h>

/// Watershed strategy selection
pso_watershed pso_watershed_select(const char *strategy)
{
    if (strcmp(strategy, "worst_so_far") == 0)
    {
        return pso_watershed_worst_so_far;
    }
    else if (strcmp(strategy, "worst_last_iter") == 0)
    {
        return pso_watershed_worst_last_iter;
    }
    else if (strcmp(strategy, "best_worst") == 0)
    {
        return pso_watershed_best_worst;
    }
    else if (strcmp(strategy, "none") == 0)
    {
        return NULL;
    }
    ZF_LOGE("Unknown watershed strategy: '%s'", strategy);
    return NULL;
}

/// Watershed strategy: minimum fitness is the worst so far
void pso_watershed_worst_so_far(PSO *pso)
{
    if (pso->worst_fitness < pso->watershed_min_fitness)
        pso->watershed_min_fitness = pso->worst_fitness;
}

// Watershed strategy: minimum fitness is the worst in last iteration
void pso_watershed_worst_last_iter(PSO *pso)
{
    pso->watershed_min_fitness = pso->worst_fitness;
}

// Watershed strategy: minimum fitness is the best of all worsts so far
void pso_watershed_best_worst(PSO *pso)
{
    if (pso->worst_fitness > pso->watershed_min_fitness)
        pso->watershed_min_fitness = pso->worst_fitness;
}
