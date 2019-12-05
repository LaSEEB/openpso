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
#include <float.h>

/// Watershed strategy selection
void pso_watershed_select(PSO *pso, const char *strategy)
{
    pso->watershed = NULL;
    if (strcmp(strategy, "worst_so_far") == 0)
    {
        // Set watershed strategy to worst so far
        pso->watershed = pso_watershed_worst_so_far;
        // Initialize watershed maximum fitness to minimum possible
        // worst so far strategy will immediately update it to a higher maximum
        pso->watershed_max_fitness = DBL_MIN;
    }
    else if (strcmp(strategy, "worst_last_iter") == 0)
    {
        // Set watershed strategy to worst in last iteration
        pso->watershed = pso_watershed_worst_last_iter;
        // No need to define initial value for watershed maximum fitness, since
        // it will always be the worst in last iteration
    }
    else if (strcmp(strategy, "best_worst") == 0)
    {
        // Set watershed strategy to best_worst
        pso->watershed = pso_watershed_best_worst;
        // Initialize watershed maximum fitness to maximum possible
        // best_worst strategy will immediately update it to a lower maximum
        pso->watershed_max_fitness = DBL_MAX;
    }
    else if (strcmp(strategy, "none") == 0)
    {
        // No watershed strategy
        pso->watershed = NULL;
    }
    else
    {
        // Unknown watershed strategy, log an error
        ZF_LOGE("Unknown watershed strategy: '%s'", strategy);
        pso->watershed = NULL;
    }
}

/// Watershed strategy: maximum fitness is the worst so far
void pso_watershed_worst_so_far(PSO *pso)
{
    if (pso->worst_fitness > pso->watershed_max_fitness)
        pso->watershed_max_fitness = pso->worst_fitness;
}

// Watershed strategy: maximum fitness is the worst in last iteration
void pso_watershed_worst_last_iter(PSO *pso)
{
    pso->watershed_max_fitness = pso->worst_fitness;
}

// Watershed strategy: maximum fitness is the best of all worsts so far
void pso_watershed_best_worst(PSO *pso)
{
    if (pso->worst_fitness < pso->watershed_max_fitness)
        pso->watershed_max_fitness = pso->worst_fitness;
}
