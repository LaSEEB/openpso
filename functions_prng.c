#include <math.h>
#include <inttypes.h>

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0



static double MRG32k3a(double state) {

	/* Separate state into sub-states. This is done in an arbitrary fashion,
	 * but must follow these rules:
	 * - s10, s11, s12 must be integers in [0, m1 - 1] and not all 0.
	 * - s20, s21, s22 must be integers in [0, m2 - 1] and not all 0. */
	double s10 = (double) ((long) (fabs(log(state)) + 1.0) % (long) m1);
	double s11 = (double) ((long) (fabs(exp(state)) + 500.0) % (long) m1);
	double s12 = (double) ((long) (sqrt(fabs(state)) + 2000.0) % (long) m1);
	double s20 = (double) ((long) (fabs(tanh(state)) + 5.0) % (long) m2);
	double s21 = (double) ((long) (fabs(log10(state)) + 750.0) % (long) m2);
	double s22 = (double) ((long) (fabs(state) + 8000.0) % (long) m2);

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

double runif01(double * vars, unsigned int nvars) {

	double seed = vars[0];
	unsigned int i = 0;
	uint64_t lseed;
	for (i = 1; i < nvars; ++i) {
		seed += vars[i] / (i + 1);
	}
	lseed = (uint64_t) seed;
	lseed ^= (lseed << 21);
	lseed ^= (lseed >> 35);
	lseed ^= (lseed << 4);
	return MRG32k3a((double) lseed * seed);
}
