/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file
 * A header library with a small hash for mixing PRNG seeds.
 *
 * @author Nuno Fachada
 */

#ifndef __MIXSEED_H_
#define __MIXSEED_H_

#include <stdint.h>

/**
 * Internal function for hashing PRNG seed with thread ID.
 *
 * @param[in] seed Global PRNG seed.
 * @param[in] tid Thread ID.
 * @return A thread-specific PRNG seed.
 */
static inline uint32_t mixseed(uint32_t seed, uint32_t tid) {

	// A different base seed for each thread
	seed += tid;

	// Hash base seed
	seed  = (seed + 0x7ed55d16) + (seed << 12);
	seed  = (seed + 0x165667b1) + (seed << 5);
	seed  = (seed ^ 0xc761c23c) ^ (seed >> 19);
	seed  = (seed + 0xd3a2646c) ^ (seed << 9);
	seed  = (seed + 0xfd7046c5) + (seed << 3);
	seed  = (seed ^ 0xb55a4f09) ^ (seed >> 16);

	// Hash thread ID
	tid = (tid ^ 61) ^ (tid >> 16);
	tid = tid + (tid << 3);
	tid = tid ^ (tid >> 4);
	tid = tid * 0x27d4eb2d;
	tid = tid ^ (tid >> 15);

	// Return mixed seed
	return seed ^ tid;
}

#endif