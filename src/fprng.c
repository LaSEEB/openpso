/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
 
#include <math.h>
#include <inttypes.h>
#include <time.h>

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

/* Separate state into sub-states. This is done in an arbitrary fashion,
 * but must follow these rules:
 * - s10, s11, s12 must be integers in [0, m1 - 1] and not all 0.
 * - s20, s21, s22 must be integers in [0, m2 - 1] and not all 0. */
static double MRG32k3a(
	double s10, double s11, double s12, double s20, double s21, double s22) {

	s10 = (double) (((long) s10) % ((long) m1));
	s11 = (double) (((long) s11) % ((long) m1));
	s12 = (double) (((long) s12) % ((long) m1));
	s20 = (double) (((long) s10) % ((long) m2));
	s21 = (double) (((long) s11) % ((long) m2));
	s22 = (double) (((long) s12) % ((long) m2));

	long k;
	double p1, p2;
	/* Component 1 */
	p1 = a12 * s11 - a13n * s10;
	k = p1 / m1;
	p1 -= k * m1;
	if (p1 < 0.0)
	p1 += m1;
	s10 = s11;
	s11 = s12;
	s12 = p1;

	/* Component 2 */
	p2 = a21 * s22 - a23n * s20;
	k = p2 / m2;
	p2 -= k * m2;
	if (p2 < 0.0)
	p2 += m2;
	s20 = s21;
	s21 = s22;
	s22 = p2;

	/* Combination */
	if (p1 <= p2)
	return ((p1 - p2 + m1) * norm);
	else
	return ((p1 - p2) * norm);
}

typedef union {
	uint64_t s64;
	struct {
		uint32_t low;
		uint32_t high;
	} s32;
	double sd;
} flex;

double runif01(double * vars, unsigned int nvars) {

	// Union variable used to mix seeds
	flex seed, pseed;

	// Maximum cycles for mixing seeds
	unsigned int maxfor = (nvars > 6) ? nvars : 6;

	// Entropy from random memory addresses
	uintptr_t masks[6] = {(uintptr_t) &seed, (uintptr_t) &MRG32k3a,
		(uintptr_t) &runif01, (uintptr_t) &vars, (uintptr_t) &nvars,
		(uintptr_t) &maxfor};

	// Seeds vector for MRG32k3a PRNG
	double seeds[6] = {masks[5], masks[3], masks[0], masks[3], masks[2], masks[1]};

	// Seed mixing
	for (unsigned int i = 0; i < maxfor; ++i) {

		// Initialize union variable with values from the given variables,
		// the iteration count and previous value in seeds vector
		seed.sd = vars[i % nvars] + i + seeds[(i + 1) % 6];

		// Add entropy from current time and memory addresses
		seed.s64 ^= ((uint64_t) time(NULL)) ^ ((uint64_t) masks[i % 6]);

		// Xor-Shift PRNG for further mixing
		seed.s64 ^= (seed.s64 << 21);
		seed.s64 ^= (seed.s64 >> 35);
		seed.s64 ^= (seed.s64 << 4);

		// Hash low bits
		seed.s32.low  = (seed.s32.low + 0x7ed55d16) + (seed.s32.low << 12);
		seed.s32.low  = (seed.s32.low + 0x165667b1) + (seed.s32.low << 5);
		seed.s32.low  = (seed.s32.low ^ 0xc761c23c) ^ (seed.s32.low >> 19);
		seed.s32.low  = (seed.s32.low + 0xd3a2646c) ^ (seed.s32.low << 9);
		seed.s32.low  = (seed.s32.low + 0xfd7046c5) + (seed.s32.low << 3);
		seed.s32.low  = (seed.s32.low ^ 0xb55a4f09) ^ (seed.s32.low >> 16);

		// Hash high bits
		seed.s32.high = (seed.s32.high ^ 61) ^ (seed.s32.high >> 16);
		seed.s32.high = seed.s32.high + (seed.s32.high << 3);
		seed.s32.high = seed.s32.high ^ (seed.s32.high >> 4);
		seed.s32.high = seed.s32.high * 0x27d4eb2d;
		seed.s32.high = seed.s32.high ^ (seed.s32.high >> 15);

		// Mix low and high bits, convert to double and add to current seed
		// vector value
		pseed.sd = seeds[i % 6];
		pseed.s32.low ^= seed.s32.high;
		pseed.s32.high ^= seed.s32.low;

		seeds[i % 6] += (double) pseed.s64;
	}

	for (unsigned int i = 0; i < 6; ++i) {
		uint32_t lseed;
		pseed.sd = (seeds[(i + 2) % 6] + seeds[(i + 3) % 6]) / 2;
		lseed = pseed.s32.high ^ pseed.s32.low;
		lseed = lseed >> (lseed >> 29);
		seeds[i] = (double) lseed;
	}

	// Return uniform value between 0 and 1 using the MRG32k3a PRNG
	return MRG32k3a(seeds[0], seeds[1], seeds[2], seeds[3], seeds[4], seeds[5]);
}
