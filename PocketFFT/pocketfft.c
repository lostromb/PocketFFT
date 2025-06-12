/*
 * This file is part of pocketfft.
 * Licensed under a 3-clause BSD style license - see LICENSE.md
 */

 /*
  *  Main implementation file.
  *
  *  Copyright (C) 2004-2018 Max-Planck-Society
  *  \author Martin Reinecke
  */

#include <math.h>
#include <string.h>
#include "NailTester.h"
#include "pocketfft.h"
#include "cmplx.h"

#define RALLOC(type,num) \
  ((type *)malloc((num)*sizeof(type)))
#define DEALLOC(ptr) \
  do { free(ptr); (ptr)=NULL; } while(0)

#define SWAP(a,b,type) \
  do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)

#ifdef __GNUC__
#define NOINLINE __attribute__((noinline))
#define WARN_UNUSED_RESULT __attribute__ ((warn_unused_result))
#else
#define NOINLINE
#define WARN_UNUSED_RESULT
#endif

const int NAIL_TEST = 0;

  // adapted from https://stackoverflow.com/questions/42792939/
  // CAUTION: this function only works for arguments in the range [-0.25; 0.25]!
static void my_sincosm1pi(double a, double* res)
{
	double s = a * a;
	/* Approximate cos(pi*x)-1 for x in [-0.25,0.25] */
	double r = -1.0369917389758117e-4;
	r = fma(r, s, 1.9294935641298806e-3);
	r = fma(r, s, -2.5806887942825395e-2);
	r = fma(r, s, 2.3533063028328211e-1);
	r = fma(r, s, -1.3352627688538006e+0);
	r = fma(r, s, 4.0587121264167623e+0);
	r = fma(r, s, -4.9348022005446790e+0);
	double c = r * s;
	/* Approximate sin(pi*x) for x in [-0.25,0.25] */
	r = 4.6151442520157035e-4;
	r = fma(r, s, -7.3700183130883555e-3);
	r = fma(r, s, 8.2145868949323936e-2);
	r = fma(r, s, -5.9926452893214921e-1);
	r = fma(r, s, 2.5501640398732688e+0);
	r = fma(r, s, -5.1677127800499516e+0);
	s = s * a;
	r = r * s;
	s = fma(a, 3.1415926535897931e+0, r);
	res[0] = c;
	res[1] = s;
}

NOINLINE static void calc_first_octant(size_t den, double* res)
{
	size_t n = (den + 4) >> 3;
	if (n == 0) return;
	res[0] = 1.; res[1] = 0.;
	if (n == 1) return;
	size_t l1 = (size_t)sqrt(n);
	for (size_t i = 1; i < l1; ++i)
		my_sincosm1pi((2. * i) / den, &res[2 * i]);
	size_t start = l1;
	while (start < n)
	{
		double cs[2];
		my_sincosm1pi((2. * start) / den, cs);
		res[2 * start] = cs[0] + 1.;
		res[2 * start + 1] = cs[1];
		size_t end = l1;
		if (start + end > n) end = n - start;
		for (size_t i = 1; i < end; ++i)
		{
			double csx[2] = { res[2 * i], res[2 * i + 1] };
			res[2 * (start + i)] = ((cs[0] * csx[0] - cs[1] * csx[1] + cs[0]) + csx[0]) + 1.;
			res[2 * (start + i) + 1] = (cs[0] * csx[1] + cs[1] * csx[0]) + cs[1] + csx[1];
		}
		start += l1;
	}
	for (size_t i = 1; i < l1; ++i)
		res[2 * i] += 1.;
}

NOINLINE static void calc_first_quadrant(size_t n, double* res)
{
	double* p = res + n;
	calc_first_octant(n << 1, p);
	size_t ndone = (n + 2) >> 2;
	size_t i = 0, idx1 = 0, idx2 = 2 * ndone - 2;
	for (; i + 1 < ndone; i += 2, idx1 += 2, idx2 -= 2)
	{
		res[idx1] = p[2 * i];
		res[idx1 + 1] = p[2 * i + 1];
		res[idx2] = p[2 * i + 3];
		res[idx2 + 1] = p[2 * i + 2];
	}
	if (i != ndone)
	{
		res[idx1] = p[2 * i];
		res[idx1 + 1] = p[2 * i + 1];
	}
}

NOINLINE static void calc_first_half(size_t n, double* res)
{
	int ndone = (n + 1) >> 1;
	double* p = res + n - 1;
	calc_first_octant(n << 2, p);
	int i4 = 0, in = n, i = 0;
	for (; i4 <= in - i4; ++i, i4 += 4) // octant 0
	{
		res[2 * i] = p[2 * i4]; res[2 * i + 1] = p[2 * i4 + 1];
	}
	for (; i4 - in <= 0; ++i, i4 += 4) // octant 1
	{
		int xm = in - i4;
		res[2 * i] = p[2 * xm + 1]; res[2 * i + 1] = p[2 * xm];
	}
	for (; i4 <= 3 * in - i4; ++i, i4 += 4) // octant 2
	{
		int xm = i4 - in;
		res[2 * i] = -p[2 * xm + 1]; res[2 * i + 1] = p[2 * xm];
	}
	for (; i < ndone; ++i, i4 += 4) // octant 3
	{
		int xm = 2 * in - i4;
		res[2 * i] = -p[2 * xm]; res[2 * i + 1] = p[2 * xm + 1];
	}
}

NOINLINE static void fill_first_quadrant(size_t n, double* res)
{
	const double hsqt2 = 0.707106781186547524400844362104849;
	size_t quart = n >> 2;
	if ((n & 7) == 0)
		res[quart] = res[quart + 1] = hsqt2;
	for (size_t i = 2, j = 2 * quart - 2; i < quart; i += 2, j -= 2)
	{
		res[j] = res[i + 1];
		res[j + 1] = res[i];
	}
}

NOINLINE static void fill_first_half(size_t n, double* res)
{
	size_t half = n >> 1;
	if ((n & 3) == 0)
		for (size_t i = 0; i < half; i += 2)
		{
			res[i + half] = -res[i + 1];
			res[i + half + 1] = res[i];
		}
	else
		for (size_t i = 2, j = 2 * half - 2; i < half; i += 2, j -= 2)
		{
			res[j] = -res[i];
			res[j + 1] = res[i + 1];
		}
}

NOINLINE static void fill_second_half(size_t n, double* res)
{
	if ((n & 1) == 0)
		for (size_t i = 0; i < n; ++i)
			res[i + n] = -res[i];
	else
		for (size_t i = 2, j = 2 * n - 2; i < n; i += 2, j -= 2)
		{
			res[j] = res[i];
			res[j + 1] = -res[i + 1];
		}
}

NOINLINE static void sincos_2pibyn_half(size_t n, double* res)
{
	if ((n & 3) == 0)
	{
		calc_first_octant(n, res);
		fill_first_quadrant(n, res);
		fill_first_half(n, res);
	}
	else if ((n & 1) == 0)
	{
		calc_first_quadrant(n, res);
		fill_first_half(n, res);
	}
	else
		calc_first_half(n, res);
}

NOINLINE static void sincos_2pibyn(size_t n, double* res)
{
	sincos_2pibyn_half(n, res);
	fill_second_half(n, res);
}

NOINLINE static size_t largest_prime_factor(size_t n)
{
	size_t res = 1;
	size_t tmp;
	while (((tmp = (n >> 1)) << 1) == n)
	{
		res = 2; n = tmp;
	}

	size_t limit = (size_t)sqrt(n + 0.01);
	for (size_t x = 3; x <= limit; x += 2)
		while (((tmp = (n / x)) * x) == n)
		{
			res = x;
			n = tmp;
			limit = (size_t)sqrt(n + 0.01);
		}
	if (n > 1) res = n;

	return res;
}

NOINLINE static double cost_guess(size_t n)
{
	const double lfp = 1.1; // penalty for non-hardcoded larger factors
	size_t ni = n;
	double result = 0.;
	size_t tmp;
	while (((tmp = (n >> 1)) << 1) == n)
	{
		result += 2; n = tmp;
	}

	size_t limit = (size_t)sqrt(n + 0.01);
	for (size_t x = 3; x <= limit; x += 2)
		while ((tmp = (n / x)) * x == n)
		{
			result += (x <= 5) ? x : lfp * x; // penalize larger prime factors
			n = tmp;
			limit = (size_t)sqrt(n + 0.01);
		}
	if (n > 1) result += (n <= 5) ? n : lfp * n;

	return result * ni;
}

/* returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n */
NOINLINE static size_t good_size(size_t n)
{
	if (n <= 6) return n;

	size_t bestfac = 2 * n;
	for (size_t f2 = 1; f2 < bestfac; f2 *= 2)
		for (size_t f23 = f2; f23 < bestfac; f23 *= 3)
			for (size_t f235 = f23; f235 < bestfac; f235 *= 5)
				for (size_t f2357 = f235; f2357 < bestfac; f2357 *= 7)
					for (size_t f235711 = f2357; f235711 < bestfac; f235711 *= 11)
						if (f235711 >= n) bestfac = f235711;
	return bestfac;
}

#define NFCT 25
typedef struct cfftp_fctdata
{
	size_t fct;
	cmplx* tw, * tws;
} cfftp_fctdata;

typedef struct cfftp_plan_i
{
	size_t length, nfct;
	cmplx* mem;
	cfftp_fctdata fct[NFCT];
} cfftp_plan_i;
typedef struct cfftp_plan_i* cfftp_plan;

#define PMC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
#define ADDC(a,b,c) { a.r=b.r+c.r; a.i=b.i+c.i; }
#define SCALEC(a,b) { a.r*=b; a.i*=b; }
#define ROT90(a) { double tmp_=a.r; a.r=-a.i; a.i=tmp_; }
#define ROTM90(a) { double tmp_=-a.r; a.r=a.i; a.i=tmp_; }
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define WA(x,i) wa[(i)-1+(x)*(ido-1)]
/* a = b*c */
#define A_EQ_B_MUL_C(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
/* a = conj(b)*c*/
#define A_EQ_CB_MUL_C(a,b,c) { a.r=b.r*c.r+b.i*c.i; a.i=b.r*c.i-b.i*c.r; }

#define PMSIGNC(a,b,c,d) { a.r=c.r+sign*d.r; a.i=c.i+sign*d.i; b.r=c.r-sign*d.r; b.i=c.i-sign*d.i; }
/* a = b*c */
#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r-sign*b.i*c.i; a.i=b.r*c.i+sign*b.i*c.r; }
/* a *= b */
#define MULPMSIGNCEQ(a,b) { double xtmp=a.r; a.r=b.r*a.r-sign*b.i*a.i; a.i=b.r*a.i+sign*b.i*xtmp; }

NOINLINE static void pass2b(size_t ido, size_t l1, const cmplx* cc, cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 2;

	if (NAIL_TEST) printf("start pass2b\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;

			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t;
				ch[(i)+ido * ((k)+l1 * (0))].r = cc[(i)+ido * ((0) + cdim * (k))].r + cc[(i)+ido * ((1) + cdim * (k))].r;
				ch[(i)+ido * ((k)+l1 * (0))].i = cc[(i)+ido * ((0) + cdim * (k))].i + cc[(i)+ido * ((1) + cdim * (k))].i;
				t.r = cc[(i)+ido * ((0) + cdim * (k))].r - cc[(i)+ido * ((1) + cdim * (k))].r;
				t.i = cc[(i)+ido * ((0) + cdim * (k))].i - cc[(i)+ido * ((1) + cdim * (k))].i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (0) * (ido - 1)].r * t.r - wa[(i)-1 + (0) * (ido - 1)].i * t.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (0) * (ido - 1)].r * t.i + wa[(i)-1 + (0) * (ido - 1)].i * t.r;
			}
		}
	}
}

NOINLINE static void pass2f(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 2;
	if (NAIL_TEST) printf("start pass2f\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
			ch[(0) + ido * ((k)+l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t;
				ch[(i)+ido * ((k)+l1 * (0))].r = cc[(i)+ido * ((0) + cdim * (k))].r + cc[(i)+ido * ((1) + cdim * (k))].r;
				ch[(i)+ido * ((k)+l1 * (0))].i = cc[(i)+ido * ((0) + cdim * (k))].i + cc[(i)+ido * ((1) + cdim * (k))].i;
				t.r = cc[(i)+ido * ((0) + cdim * (k))].r - cc[(i)+ido * ((1) + cdim * (k))].r;
				t.i = cc[(i)+ido * ((0) + cdim * (k))].i - cc[(i)+ido * ((1) + cdim * (k))].i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (0) * (ido - 1)].r * t.r + wa[(i)-1 + (0) * (ido - 1)].i * t.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (0) * (ido - 1)].r * t.i - wa[(i)-1 + (0) * (ido - 1)].i * t.r;
			}
		}
	}
}

NOINLINE static void pass3b(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 3;
	const double tw1r = -0.5, tw1i = 0.86602540378443864676;

	if (NAIL_TEST) printf("start pass3b\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
			t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i;
			cmplx ca, cb;
			ca.r = t0.r + tw1r * t1.r; ca.i = t0.i + tw1r * t1.i;
			cb.i = tw1i * t2.r;
			cb.r = -(tw1i * t2.i);
			ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = ca.r - cb.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = ca.i - cb.i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
			t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i;
			cmplx ca, cb;
			ca.r = t0.r + tw1r * t1.r;
			ca.i = t0.i + tw1r * t1.i;
			cb.i = tw1i * t2.r;
			cb.r = -(tw1i * t2.i);
			ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = ca.r - cb.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = ca.i - cb.i;

			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t0 = cc[(i)+ido * ((0) + cdim * (k))], t1, t2;
				t1.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((2) + cdim * (k))].r;
				t1.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((2) + cdim * (k))].i;
				t2.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((2) + cdim * (k))].r;
				t2.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((2) + cdim * (k))].i;
				ch[(i)+ido * ((k)+l1 * (0))].r = t0.r + t1.r;
				ch[(i)+ido * ((k)+l1 * (0))].i = t0.i + t1.i;
				cmplx ca, cb, da, db;
				ca.r = t0.r + tw1r * t1.r;
				ca.i = t0.i + tw1r * t1.i;
				cb.i = tw1i * t2.r;
				cb.r = -(tw1i * t2.i);
				da.r = ca.r + cb.r;
				da.i = ca.i + cb.i;
				db.r = ca.r - cb.r;
				db.i = ca.i - cb.i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
				ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * db.r - wa[(i)-1 + (2 - 1) * (ido - 1)].i * db.i;
				ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * db.i + wa[(i)-1 + (2 - 1) * (ido - 1)].i * db.r;
			}
		}
	}
}

NOINLINE static void pass3f(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 3;
	const double tw1r = -0.5, tw1i = -0.86602540378443864676;

	if (NAIL_TEST) printf("start pass3f\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
			t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i;
			cmplx ca, cb;
			ca.r = t0.r + tw1r * t1.r;
			ca.i = t0.i + tw1r * t1.i;
			cb.i = tw1i * t2.r;
			cb.r = -(tw1i * t2.i);
			ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = ca.r - cb.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = ca.i - cb.i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
			t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i;
			cmplx ca, cb;
			ca.r = t0.r + tw1r * t1.r;
			ca.i = t0.i + tw1r * t1.i;
			cb.i = tw1i * t2.r;
			cb.r = -(tw1i * t2.i);
			ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = ca.r - cb.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = ca.i - cb.i;

			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t3 = cc[(i)+ido * ((0) + cdim * (k))], t4, t5;
				t4.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((2) + cdim * (k))].r;
				t4.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((2) + cdim * (k))].i;
				t5.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((2) + cdim * (k))].r;
				t5.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((2) + cdim * (k))].i;
				ch[(i)+ido * ((k)+l1 * (0))].r = t3.r + t4.r;
				ch[(i)+ido * ((k)+l1 * (0))].i = t3.i + t4.i;
				cmplx ca2, cb2, da, db;
				ca2.r = t3.r + tw1r * t4.r;
				ca2.i = t3.i + tw1r * t4.i;
				cb2.i = tw1i * t5.r;
				cb2.r = -(tw1i * t5.i);
				da.r = ca2.r + cb2.r;
				da.i = ca2.i + cb2.i;
				db.r = ca2.r - cb2.r;
				db.i = ca2.i - cb2.i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
				ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * db.r + wa[(i)-1 + (2 - 1) * (ido - 1)].i * db.i;
				ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * db.i - wa[(i)-1 + (2 - 1) * (ido - 1)].i * db.r;
			}
		}
	}
}

NOINLINE static void pass4b(size_t ido, size_t l1, const cmplx* cc, cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 4;
	if (NAIL_TEST) printf("start pass4b\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1, t2, t3, t4;
			t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
			t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
			t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
			t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			ROT90(t4)
			ch[(0) + ido * ((k)+l1 * (0))].r = t2.r + t3.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t2.i + t3.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = t2.r - t3.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = t2.i - t3.i;
			ch[(0) + ido * ((k)+l1 * (1))].r = t1.r + t4.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = t1.i + t4.i;
			ch[(0) + ido * ((k)+l1 * (3))].r = t1.r - t4.r;
			ch[(0) + ido * ((k)+l1 * (3))].i = t1.i - t4.i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1, t2, t3, t4;
			t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
			t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
			t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
			t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			ROT90(t4)
				ch[(0) + ido * ((k)+l1 * (0))].r = t2.r + t3.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t2.i + t3.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = t2.r - t3.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = t2.i - t3.i;
			ch[(0) + ido * ((k)+l1 * (1))].r = t1.r + t4.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = t1.i + t4.i;
			ch[(0) + ido * ((k)+l1 * (3))].r = t1.r - t4.r;
			ch[(0) + ido * ((k)+l1 * (3))].i = t1.i - t4.i;

			for (size_t i = 1; i < ido; ++i)
			{
				cmplx c2, c3, c4, t5, t6, t7, t8;
				cmplx cc0 = cc[(i)+ido * ((0) + cdim * (k))],
					cc1 = cc[(i)+ido * ((1) + cdim * (k))],
					cc2 = cc[(i)+ido * ((2) + cdim * (k))],
					cc3 = cc[(i)+ido * ((3) + cdim * (k))];
				t6.r = cc0.r + cc2.r;
				t6.i = cc0.i + cc2.i;
				t5.r = cc0.r - cc2.r;
				t5.i = cc0.i - cc2.i;
				t7.r = cc1.r + cc3.r;
				t7.i = cc1.i + cc3.i;
				t8.r = cc1.r - cc3.r;
				t8.i = cc1.i - cc3.i;
				ROT90(t8)
				cmplx wa0 = wa[(i)-1 + (0) * (ido - 1)],
					wa1 = wa[(i)-1 + (1) * (ido - 1)],
					wa2 = wa[(i)-1 + (2) * (ido - 1)];
				ch[(i)+ido * ((k)+l1 * (0))].r = t6.r + t7.r;
				ch[(i)+ido * ((k)+l1 * (0))].i = t6.i + t7.i;
				c3.r = t6.r - t7.r;
				c3.i = t6.i - t7.i;
				c2.r = t5.r + t8.r;
				c2.i = t5.i + t8.i;
				c4.r = t5.r - t8.r;
				c4.i = t5.i - t8.i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa0.r * c2.r - wa0.i * c2.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa0.r * c2.i + wa0.i * c2.r;
				ch[(i)+ido * ((k)+l1 * (2))].r = wa1.r * c3.r - wa1.i * c3.i;
				ch[(i)+ido * ((k)+l1 * (2))].i = wa1.r * c3.i + wa1.i * c3.r;
				ch[(i)+ido * ((k)+l1 * (3))].r = wa2.r * c4.r - wa2.i * c4.i;
				ch[(i)+ido * ((k)+l1 * (3))].i = wa2.r * c4.i + wa2.i * c4.r;
			}
		}
	}
}
NOINLINE static void pass4f(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 4;

	if (NAIL_TEST) printf("start pass4f\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1, t2, t3, t4;
			t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
			t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
			t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
			t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			ROTM90(t4)
			ch[(0) + ido * ((k)+l1 * (0))].r = t2.r + t3.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t2.i + t3.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = t2.r - t3.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = t2.i - t3.i;
			ch[(0) + ido * ((k)+l1 * (1))].r = t1.r + t4.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = t1.i + t4.i;
			ch[(0) + ido * ((k)+l1 * (3))].r = t1.r - t4.r;
			ch[(0) + ido * ((k)+l1 * (3))].i = t1.i - t4.i;
		}
	}
	else
	{
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1, t2, t3, t4;
			t2.r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
			t2.i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
			t1.r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
			t1.i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
			t3.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r;
			t3.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i;
			t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r;
			t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			ROTM90(t4)
			ch[(0) + ido * ((k)+l1 * (0))].r = t2.r + t3.r;
			ch[(0) + ido * ((k)+l1 * (0))].i = t2.i + t3.i;
			ch[(0) + ido * ((k)+l1 * (2))].r = t2.r - t3.r;
			ch[(0) + ido * ((k)+l1 * (2))].i = t2.i - t3.i;
			ch[(0) + ido * ((k)+l1 * (1))].r = t1.r + t4.r;
			ch[(0) + ido * ((k)+l1 * (1))].i = t1.i + t4.i;
			ch[(0) + ido * ((k)+l1 * (3))].r = t1.r - t4.r;
			ch[(0) + ido * ((k)+l1 * (3))].i = t1.i - t4.i;
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx c2, c3, c4, t5, t6, t7, t8;
				cmplx cc0 = cc[(i)+ido * ((0) + cdim * (k))],
					cc1 = cc[(i)+ido * ((1) + cdim * (k))],
					cc2 = cc[(i)+ido * ((2) + cdim * (k))],
					cc3 = cc[(i)+ido * ((3) + cdim * (k))];
				t6.r = cc0.r + cc2.r;
				t6.i = cc0.i + cc2.i;
				t5.r = cc0.r - cc2.r;
				t5.i = cc0.i - cc2.i;
				t7.r = cc1.r + cc3.r;
				t7.i = cc1.i + cc3.i;
				t8.r = cc1.r - cc3.r;
				t8.i = cc1.i - cc3.i;
				ROTM90(t8)
				cmplx wa0 = wa[(i)-1 + (0) * (ido - 1)],
					wa1 = wa[(i)-1 + (1) * (ido - 1)],
					wa2 = wa[(i)-1 + (2) * (ido - 1)];
				ch[(i)+ido * ((k)+l1 * (0))].r = t6.r + t7.r;
				ch[(i)+ido * ((k)+l1 * (0))].i = t6.i + t7.i;
				c3.r = t6.r - t7.r;
				c3.i = t6.i - t7.i;
				c2.r = t5.r + t8.r;
				c2.i = t5.i + t8.i;
				c4.r = t5.r - t8.r; 
				c4.i = t5.i - t8.i;
				ch[(i)+ido * ((k)+l1 * (1))].r = wa0.r * c2.r + wa0.i * c2.i;
				ch[(i)+ido * ((k)+l1 * (1))].i = wa0.r * c2.i - wa0.i * c2.r;
				ch[(i)+ido * ((k)+l1 * (2))].r = wa1.r * c3.r + wa1.i * c3.i;
				ch[(i)+ido * ((k)+l1 * (2))].i = wa1.r * c3.i - wa1.i * c3.r;
				ch[(i)+ido * ((k)+l1 * (3))].r = wa2.r * c4.r + wa2.i * c4.i;
				ch[(i)+ido * ((k)+l1 * (3))].i = wa2.r * c4.i - wa2.i * c4.r;
			}
		}
	}
}

#define PREP5(idx) \
        cmplx t0 = CC(idx,0,k), t1, t2, t3, t4; \
        PMC (t1,t4,CC(idx,1,k),CC(idx,4,k)) \
        PMC (t2,t3,CC(idx,2,k),CC(idx,3,k)) \
        CH(idx,k,0).r=t0.r+t1.r+t2.r; \
        CH(idx,k,0).i=t0.i+t1.i+t2.i;

#define PARTSTEP5a(u1,u2,twar,twbr,twai,twbi) \
        { \
        cmplx ca,cb; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(CH(0,k,u1),CH(0,k,u2),ca,cb) \
        }

#define PARTSTEP5b(u1,u2,twar,twbr,twai,twbi) \
        { \
        cmplx ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb) \
        A_EQ_B_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_B_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
NOINLINE static void pass5b(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 5;
	const double tw1r = 0.3090169943749474241,
		tw1i = 0.95105651629515357212,
		tw2r = -0.8090169943749474241,
		tw2i = 0.58778525229247312917;

	if (NAIL_TEST) printf("start pass5b\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
				t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
			} {
				t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			} ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
			{
				cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
					ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
					ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i - cb.i;
				}
			}
		}
	else
		for (size_t k = 0; k < l1; ++k)
		{
			{
				cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
					t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
				} {
					t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
				} ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
				{
					cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
						ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
						ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i - cb.i;
					}
				}
			}
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t0 = cc[(i)+ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
					t1.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((4) + cdim * (k))].r; t1.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((4) + cdim * (k))].i; t4.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((4) + cdim * (k))].r; t4.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((4) + cdim * (k))].i;
				} {
					t2.r = cc[(i)+ido * ((2) + cdim * (k))].r + cc[(i)+ido * ((3) + cdim * (k))].r; t2.i = cc[(i)+ido * ((2) + cdim * (k))].i + cc[(i)+ido * ((3) + cdim * (k))].i; t3.r = cc[(i)+ido * ((2) + cdim * (k))].r - cc[(i)+ido * ((3) + cdim * (k))].r; t3.i = cc[(i)+ido * ((2) + cdim * (k))].i - cc[(i)+ido * ((3) + cdim * (k))].i;
				} ch[(i)+ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i)+ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
				{
					cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
						da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
					} {
						ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (4))].r = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.r - wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (4))].i = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.i + wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
						da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
					} {
						ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.r - wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.i + wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (3))].r = wa[(i)-1 + (3 - 1) * (ido - 1)].r * db.r - wa[(i)-1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (3))].i = wa[(i)-1 + (3 - 1) * (ido - 1)].r * db.i + wa[(i)-1 + (3 - 1) * (ido - 1)].i * db.r;
					}
				}
			}
		}
}
#define PARTSTEP5f(u1,u2,twar,twbr,twai,twbi) \
        { \
        cmplx ca,cb,da,db; \
        ca.r=t0.r+twar*t1.r+twbr*t2.r; \
        ca.i=t0.i+twar*t1.i+twbr*t2.i; \
        cb.i=twai*t4.r twbi*t3.r; \
        cb.r=-(twai*t4.i twbi*t3.i); \
        PMC(da,db,ca,cb) \
        A_EQ_CB_MUL_C (CH(i,k,u1),WA(u1-1,i),da) \
        A_EQ_CB_MUL_C (CH(i,k,u2),WA(u2-1,i),db) \
        }
NOINLINE static void pass5f(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa)
{
	const size_t cdim = 5;
	const double tw1r = 0.3090169943749474241,
		tw1i = -0.95105651629515357212,
		tw2r = -0.8090169943749474241,
		tw2i = -0.58778525229247312917;

	if (NAIL_TEST) printf("start pass5f\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
				t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
			} {
				t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
			} ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
			{
				cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
					ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
				}
			}
				{
					cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
						ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i - cb.i;
					}
				}
		}
	else
		for (size_t k = 0; k < l1; ++k)
		{
			{
				cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
					t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
				} {
					t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
				} ch[(0) + ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
				{
					cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
						ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
						ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i - cb.i;
					}
				}
			}
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t0 = cc[(i)+ido * ((0) + cdim * (k))], t1, t2, t3, t4; {
					t1.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((4) + cdim * (k))].r; t1.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((4) + cdim * (k))].i; t4.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((4) + cdim * (k))].r; t4.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((4) + cdim * (k))].i;
				} {
					t2.r = cc[(i)+ido * ((2) + cdim * (k))].r + cc[(i)+ido * ((3) + cdim * (k))].r; t2.i = cc[(i)+ido * ((2) + cdim * (k))].i + cc[(i)+ido * ((3) + cdim * (k))].i; t3.r = cc[(i)+ido * ((2) + cdim * (k))].r - cc[(i)+ido * ((3) + cdim * (k))].r; t3.i = cc[(i)+ido * ((2) + cdim * (k))].i - cc[(i)+ido * ((3) + cdim * (k))].i;
				} ch[(i)+ido * ((k)+l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i)+ido * ((k)+l1 * (0))].i = t0.i + t1.i + t2.i;
				{
					cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i); {
						da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
					} {
						ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (4))].r = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.r + wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (4))].i = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.i - wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.r;
					}
				}
					{
						cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						} {
							ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.r + wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.i - wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.r;
						} {
							ch[(i)+ido * ((k)+l1 * (3))].r = wa[(i)-1 + (3 - 1) * (ido - 1)].r * db.r + wa[(i)-1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (3))].i = wa[(i)-1 + (3 - 1) * (ido - 1)].r * db.i - wa[(i)-1 + (3 - 1) * (ido - 1)].i * db.r;
						}
					}
			}
		}
}

#define PREP7(idx) \
        cmplx t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7; \
        PMC (t2,t7,CC(idx,1,k),CC(idx,6,k)) \
        PMC (t3,t6,CC(idx,2,k),CC(idx,5,k)) \
        PMC (t4,t5,CC(idx,3,k),CC(idx,4,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i;

#define PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,out1,out2) \
        { \
        cmplx ca,cb; \
        ca.r=t1.r+x1*t2.r+x2*t3.r+x3*t4.r; \
        ca.i=t1.i+x1*t2.i+x2*t3.i+x3*t4.i; \
        cb.i=y1*t7.r y2*t6.r y3*t5.r; \
        cb.r=-(y1*t7.i y2*t6.i y3*t5.i); \
        PMC(out1,out2,ca,cb) \
        }
#define PARTSTEP7a(u1,u2,x1,x2,x3,y1,y2,y3) \
        PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,CH(0,k,u1),CH(0,k,u2))
#define PARTSTEP7(u1,u2,x1,x2,x3,y1,y2,y3) \
        { \
        cmplx da,db; \
        PARTSTEP7a0(u1,u2,x1,x2,x3,y1,y2,y3,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void pass7(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa, const int sign)
{
	const size_t cdim = 7;
	const double tw1r = 0.623489801858733530525,
		tw1i = sign * 0.7818314824680298087084,
		tw2r = -0.222520933956314404289,
		tw2i = sign * 0.9749279121818236070181,
		tw3r = -0.9009688679024191262361,
		tw3i = sign * 0.4338837391175581204758;

	if (NAIL_TEST) printf("start pass7\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7; {
				t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
			} {
				t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
			} {
				t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
			} ch[(0) + ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
			{
				cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i); {
					ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (6))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i); {
					ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (5))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i); {
					ch[(0) + ido * ((k)+l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
				}
			}
		}
	else
		for (size_t k = 0; k < l1; ++k)
		{
			{
				cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7; {
					t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
				} {
					t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
				} {
					t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
				} ch[(0) + ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
				{
					cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i); {
						ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (6))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i); {
						ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (5))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i); {
						ch[(0) + ido * ((k)+l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i - cb.i;
					}
				}
			}
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t1 = cc[(i)+ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7; {
					t2.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((6) + cdim * (k))].r; t2.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((6) + cdim * (k))].i; t7.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((6) + cdim * (k))].r; t7.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((6) + cdim * (k))].i;
				} {
					t3.r = cc[(i)+ido * ((2) + cdim * (k))].r + cc[(i)+ido * ((5) + cdim * (k))].r; t3.i = cc[(i)+ido * ((2) + cdim * (k))].i + cc[(i)+ido * ((5) + cdim * (k))].i; t6.r = cc[(i)+ido * ((2) + cdim * (k))].r - cc[(i)+ido * ((5) + cdim * (k))].r; t6.i = cc[(i)+ido * ((2) + cdim * (k))].i - cc[(i)+ido * ((5) + cdim * (k))].i;
				} {
					t4.r = cc[(i)+ido * ((3) + cdim * (k))].r + cc[(i)+ido * ((4) + cdim * (k))].r; t4.i = cc[(i)+ido * ((3) + cdim * (k))].i + cc[(i)+ido * ((4) + cdim * (k))].i; t5.r = cc[(i)+ido * ((3) + cdim * (k))].r - cc[(i)+ido * ((4) + cdim * (k))].r; t5.i = cc[(i)+ido * ((3) + cdim * (k))].i - cc[(i)+ido * ((4) + cdim * (k))].i;
				} ch[(i)+ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(i)+ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (6))].r = wa[(i)-1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (6))].i = wa[(i)-1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (6 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (5))].r = wa[(i)-1 + (5 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (5 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (5))].i = wa[(i)-1 + (5 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (5 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (3))].r = wa[(i)-1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (3))].i = wa[(i)-1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (3 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (4))].r = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (4))].i = wa[(i)-1 + (4 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (4 - 1) * (ido - 1)].i * db.r;
					}
				}
			}
		}
}

#define PREP11(idx) \
        cmplx t1 = CC(idx,0,k), t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; \
        PMC (t2,t11,CC(idx,1,k),CC(idx,10,k)) \
        PMC (t3,t10,CC(idx,2,k),CC(idx, 9,k)) \
        PMC (t4,t9 ,CC(idx,3,k),CC(idx, 8,k)) \
        PMC (t5,t8 ,CC(idx,4,k),CC(idx, 7,k)) \
        PMC (t6,t7 ,CC(idx,5,k),CC(idx, 6,k)) \
        CH(idx,k,0).r=t1.r+t2.r+t3.r+t4.r+t5.r+t6.r; \
        CH(idx,k,0).i=t1.i+t2.i+t3.i+t4.i+t5.i+t6.i;

#define PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,out1,out2) \
        { \
        cmplx ca,cb; \
        ca.r=t1.r+x1*t2.r+x2*t3.r+x3*t4.r+x4*t5.r+x5*t6.r; \
        ca.i=t1.i+x1*t2.i+x2*t3.i+x3*t4.i+x4*t5.i+x5*t6.i; \
        cb.i=y1*t11.r y2*t10.r y3*t9.r y4*t8.r y5*t7.r; \
        cb.r=-(y1*t11.i y2*t10.i y3*t9.i y4*t8.i y5*t7.i ); \
        PMC(out1,out2,ca,cb) \
        }
#define PARTSTEP11a(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,CH(0,k,u1),CH(0,k,u2))
#define PARTSTEP11(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) \
        { \
        cmplx da,db; \
        PARTSTEP11a0(u1,u2,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,da,db) \
        MULPMSIGNC (CH(i,k,u1),WA(u1-1,i),da) \
        MULPMSIGNC (CH(i,k,u2),WA(u2-1,i),db) \
        }

NOINLINE static void pass11(size_t ido, size_t l1, const cmplx* cc,
	cmplx* ch, const cmplx* wa, const int sign)
{
	const size_t cdim = 11;
	const double tw1r = 0.8412535328311811688618,
		tw1i = sign * 0.5406408174555975821076,
		tw2r = 0.4154150130018864255293,
		tw2i = sign * 0.9096319953545183714117,
		tw3r = -0.1423148382732851404438,
		tw3i = sign * 0.9898214418809327323761,
		tw4r = -0.6548607339452850640569,
		tw4i = sign * 0.755749574354258283774,
		tw5r = -0.9594929736144973898904,
		tw5i = sign * 0.2817325568414296977114;

	if (NAIL_TEST) printf("start pass11\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	if (ido == 1)
		for (size_t k = 0; k < l1; ++k)
		{
			cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; {
				t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
			} {
				t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
			} {
				t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
			} {
				t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
			} {
				t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
			} ch[(0) + ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
			{
				cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i); {
					ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (10))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i); {
					ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (9))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i); {
					ch[(0) + ido * ((k)+l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (8))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i); {
					ch[(0) + ido * ((k)+l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (7))].i = ca.i - cb.i;
				}
			}
			{
				cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i); {
					ch[(0) + ido * ((k)+l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (6))].i = ca.i - cb.i;
				}
			}
		}
	else
		for (size_t k = 0; k < l1; ++k)
		{
			{
				cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; {
					t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
				} {
					t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
				} {
					t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
				} {
					t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
				} {
					t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
				} ch[(0) + ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
				{
					cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i); {
						ch[(0) + ido * ((k)+l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (10))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i); {
						ch[(0) + ido * ((k)+l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (9))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i); {
						ch[(0) + ido * ((k)+l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (8))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i); {
						ch[(0) + ido * ((k)+l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (7))].i = ca.i - cb.i;
					}
				}
				{
					cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i); {
						ch[(0) + ido * ((k)+l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k)+l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k)+l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k)+l1 * (6))].i = ca.i - cb.i;
					}
				}
			}
			for (size_t i = 1; i < ido; ++i)
			{
				cmplx t1 = cc[(i)+ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11; {
					t2.r = cc[(i)+ido * ((1) + cdim * (k))].r + cc[(i)+ido * ((10) + cdim * (k))].r; t2.i = cc[(i)+ido * ((1) + cdim * (k))].i + cc[(i)+ido * ((10) + cdim * (k))].i; t11.r = cc[(i)+ido * ((1) + cdim * (k))].r - cc[(i)+ido * ((10) + cdim * (k))].r; t11.i = cc[(i)+ido * ((1) + cdim * (k))].i - cc[(i)+ido * ((10) + cdim * (k))].i;
				} {
					t3.r = cc[(i)+ido * ((2) + cdim * (k))].r + cc[(i)+ido * ((9) + cdim * (k))].r; t3.i = cc[(i)+ido * ((2) + cdim * (k))].i + cc[(i)+ido * ((9) + cdim * (k))].i; t10.r = cc[(i)+ido * ((2) + cdim * (k))].r - cc[(i)+ido * ((9) + cdim * (k))].r; t10.i = cc[(i)+ido * ((2) + cdim * (k))].i - cc[(i)+ido * ((9) + cdim * (k))].i;
				} {
					t4.r = cc[(i)+ido * ((3) + cdim * (k))].r + cc[(i)+ido * ((8) + cdim * (k))].r; t4.i = cc[(i)+ido * ((3) + cdim * (k))].i + cc[(i)+ido * ((8) + cdim * (k))].i; t9.r = cc[(i)+ido * ((3) + cdim * (k))].r - cc[(i)+ido * ((8) + cdim * (k))].r; t9.i = cc[(i)+ido * ((3) + cdim * (k))].i - cc[(i)+ido * ((8) + cdim * (k))].i;
				} {
					t5.r = cc[(i)+ido * ((4) + cdim * (k))].r + cc[(i)+ido * ((7) + cdim * (k))].r; t5.i = cc[(i)+ido * ((4) + cdim * (k))].i + cc[(i)+ido * ((7) + cdim * (k))].i; t8.r = cc[(i)+ido * ((4) + cdim * (k))].r - cc[(i)+ido * ((7) + cdim * (k))].r; t8.i = cc[(i)+ido * ((4) + cdim * (k))].i - cc[(i)+ido * ((7) + cdim * (k))].i;
				} {
					t6.r = cc[(i)+ido * ((5) + cdim * (k))].r + cc[(i)+ido * ((6) + cdim * (k))].r; t6.i = cc[(i)+ido * ((5) + cdim * (k))].i + cc[(i)+ido * ((6) + cdim * (k))].i; t7.r = cc[(i)+ido * ((5) + cdim * (k))].r - cc[(i)+ido * ((6) + cdim * (k))].r; t7.i = cc[(i)+ido * ((5) + cdim * (k))].i - cc[(i)+ido * ((6) + cdim * (k))].i;
				} ch[(i)+ido * ((k)+l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(i)+ido * ((k)+l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (1))].r = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (1))].i = wa[(i)-1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (1 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (10))].r = wa[(i)-1 + (10 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (10 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (10))].i = wa[(i)-1 + (10 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (10 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (2))].r = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (2))].i = wa[(i)-1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (2 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (9))].r = wa[(i)-1 + (9 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (9 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (9))].i = wa[(i)-1 + (9 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (9 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (3))].r = wa[(i)-1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (3))].i = wa[(i)-1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (3 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (8))].r = wa[(i)-1 + (8 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (8 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (8))].i = wa[(i)-1 + (8 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (8 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (4))].r = wa[(i)-1 + (4 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (4 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (4))].i = wa[(i)-1 + (4 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (4 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (7))].r = wa[(i)-1 + (7 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (7 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (7))].i = wa[(i)-1 + (7 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (7 - 1) * (ido - 1)].i * db.r;
					}
				}
				{
					cmplx da, db; {
						cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i); {
							da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
						}
					} {
						ch[(i)+ido * ((k)+l1 * (5))].r = wa[(i)-1 + (5 - 1) * (ido - 1)].r * da.r - sign * wa[(i)-1 + (5 - 1) * (ido - 1)].i * da.i; ch[(i)+ido * ((k)+l1 * (5))].i = wa[(i)-1 + (5 - 1) * (ido - 1)].r * da.i + sign * wa[(i)-1 + (5 - 1) * (ido - 1)].i * da.r;
					} {
						ch[(i)+ido * ((k)+l1 * (6))].r = wa[(i)-1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i)-1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i)+ido * ((k)+l1 * (6))].i = wa[(i)-1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i)-1 + (6 - 1) * (ido - 1)].i * db.r;
					}
				}
			}
		}
}

#define CX(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CX2(a,b) cc[(a)+idl1*(b)]
#define CH2(a,b) ch[(a)+idl1*(b)]

NOINLINE static int passg(size_t ido, size_t ip, size_t l1,
	cmplx* cc, cmplx* ch, const cmplx* wa,
	const cmplx* csarr, const int sign)
{
	const size_t cdim = ip;
	size_t ipph = (ip + 1) / 2;
	size_t idl1 = ido * l1;

	if (NAIL_TEST) printf("start passg\n");
	if (NAIL_TEST) NailTestPrintComplexArray(cc, cdim * l1);
	cmplx* wal = RALLOC(cmplx, ip);
	if (!wal) return -1;
	wal[0] = (cmplx){ 1.,0. };		
	for (size_t i = 1; i < ip; ++i)
		wal[i] = (cmplx){ csarr[i].r,sign * csarr[i].i };

	for (size_t k = 0; k < l1; ++k)
		for (size_t i = 0; i < ido; ++i)
			ch[(i)+ido * ((k)+l1 * (0))] = cc[(i)+ido * ((0) + cdim * (k))];
	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)
		for (size_t k = 0; k < l1; ++k)
			for (size_t i = 0; i < ido; ++i)
			{
				ch[(i)+ido * ((k)+l1 * (j))].r = cc[(i)+ido * ((j)+cdim * (k))].r + cc[(i)+ido * ((jc)+cdim * (k))].r; ch[(i)+ido * ((k)+l1 * (j))].i = cc[(i)+ido * ((j)+cdim * (k))].i + cc[(i)+ido * ((jc)+cdim * (k))].i; ch[(i)+ido * ((k)+l1 * (jc))].r = cc[(i)+ido * ((j)+cdim * (k))].r - cc[(i)+ido * ((jc)+cdim * (k))].r; ch[(i)+ido * ((k)+l1 * (jc))].i = cc[(i)+ido * ((j)+cdim * (k))].i - cc[(i)+ido * ((jc)+cdim * (k))].i;
			}
				for (size_t k = 0; k < l1; ++k)
					for (size_t i = 0; i < ido; ++i)
					{
						cmplx tmp = ch[(i)+ido * ((k)+l1 * (0))];
						for (size_t j = 1; j < ipph; ++j)
						{
							tmp.r = tmp.r + ch[(i)+ido * ((k)+l1 * (j))].r; tmp.i = tmp.i + ch[(i)+ido * ((k)+l1 * (j))].i;
						}
						cc[(i)+ido * ((k)+l1 * (0))] = tmp;
					}
	for (size_t l = 1, lc = ip - 1; l < ipph; ++l, --lc)
	{
		// j=0
		for (size_t ik = 0; ik < idl1; ++ik)
		{
			cc[(ik)+idl1 * (l)].r = ch[(ik)+idl1 * (0)].r + wal[l].r * ch[(ik)+idl1 * (1)].r + wal[2 * l].r * ch[(ik)+idl1 * (2)].r;
			cc[(ik)+idl1 * (l)].i = ch[(ik)+idl1 * (0)].i + wal[l].r * ch[(ik)+idl1 * (1)].i + wal[2 * l].r * ch[(ik)+idl1 * (2)].i;
			cc[(ik)+idl1 * (lc)].r = -wal[l].i * ch[(ik)+idl1 * (ip - 1)].i - wal[2 * l].i * ch[(ik)+idl1 * (ip - 2)].i;
			cc[(ik)+idl1 * (lc)].i = wal[l].i * ch[(ik)+idl1 * (ip - 1)].r + wal[2 * l].i * ch[(ik)+idl1 * (ip - 2)].r;
		}

		size_t iwal = 2 * l;
		size_t j = 3, jc = ip - 3;
		for (; j < ipph - 1; j += 2, jc -= 2)
		{
			iwal += l; if (iwal > ip) iwal -= ip;
			cmplx xwal = wal[iwal];
			iwal += l; if (iwal > ip) iwal -= ip;
			cmplx xwal2 = wal[iwal];
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cc[(ik)+idl1 * (l)].r += ch[(ik)+idl1 * (j)].r * xwal.r + ch[(ik)+idl1 * (j + 1)].r * xwal2.r;
				cc[(ik)+idl1 * (l)].i += ch[(ik)+idl1 * (j)].i * xwal.r + ch[(ik)+idl1 * (j + 1)].i * xwal2.r;
				cc[(ik)+idl1 * (lc)].r -= ch[(ik)+idl1 * (jc)].i * xwal.i + ch[(ik)+idl1 * (jc - 1)].i * xwal2.i;
				cc[(ik)+idl1 * (lc)].i += ch[(ik)+idl1 * (jc)].r * xwal.i + ch[(ik)+idl1 * (jc - 1)].r * xwal2.i;
			}
		}
		for (; j < ipph; ++j, --jc)
		{
			iwal += l; if (iwal > ip) iwal -= ip;
			cmplx xwal = wal[iwal];
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cc[(ik)+idl1 * (l)].r += ch[(ik)+idl1 * (j)].r * xwal.r;
				cc[(ik)+idl1 * (l)].i += ch[(ik)+idl1 * (j)].i * xwal.r;
				cc[(ik)+idl1 * (lc)].r -= ch[(ik)+idl1 * (jc)].i * xwal.i;
				cc[(ik)+idl1 * (lc)].i += ch[(ik)+idl1 * (jc)].r * xwal.i;
			}
		}
	}
	DEALLOC(wal);

	// shuffling and twiddling
	if (ido == 1)
		for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cmplx t1 = cc[(ik)+idl1 * (j)], t2 = cc[(ik)+idl1 * (jc)];
				{ cc[(ik)+idl1 * (j)].r = t1.r + t2.r; cc[(ik)+idl1 * (j)].i = t1.i + t2.i; cc[(ik)+idl1 * (jc)].r = t1.r - t2.r; cc[(ik)+idl1 * (jc)].i = t1.i - t2.i; }
			}
	else
	{
		for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)
			for (size_t k = 0; k < l1; ++k)
			{
				cmplx t1 = cc[(0) + ido * ((k)+l1 * (j))], t2 = cc[(0) + ido * ((k)+l1 * (jc))];
				{ cc[(0) + ido * ((k)+l1 * (j))].r = t1.r + t2.r; cc[(0) + ido * ((k)+l1 * (j))].i = t1.i + t2.i; cc[(0) + ido * ((k)+l1 * (jc))].r = t1.r - t2.r; cc[(0) + ido * ((k)+l1 * (jc))].i = t1.i - t2.i; }
					for (size_t i = 1; i < ido; ++i)
					{
						cmplx x1, x2;
						{ x1.r = cc[(i)+ido * ((k)+l1 * (j))].r + cc[(i)+ido * ((k)+l1 * (jc))].r; x1.i = cc[(i)+ido * ((k)+l1 * (j))].i + cc[(i)+ido * ((k)+l1 * (jc))].i; x2.r = cc[(i)+ido * ((k)+l1 * (j))].r - cc[(i)+ido * ((k)+l1 * (jc))].r; x2.i = cc[(i)+ido * ((k)+l1 * (j))].i - cc[(i)+ido * ((k)+l1 * (jc))].i; }
							size_t idij = (j - 1) * (ido - 1) + i - 1;
							{ cc[(i)+ido * ((k)+l1 * (j))].r = wa[idij].r * x1.r - sign * wa[idij].i * x1.i; cc[(i)+ido * ((k)+l1 * (j))].i = wa[idij].r * x1.i + sign * wa[idij].i * x1.r; }
							idij = (jc - 1) * (ido - 1) + i - 1;
							{ cc[(i)+ido * ((k)+l1 * (jc))].r = wa[idij].r * x2.r - sign * wa[idij].i * x2.i; cc[(i)+ido * ((k)+l1 * (jc))].i = wa[idij].r * x2.i + sign * wa[idij].i * x2.r; }
					}
			}
	}
	return 0;
}

#undef CH2
#undef CX2
#undef CX

NOINLINE WARN_UNUSED_RESULT static int pass_all(cfftp_plan plan, cmplx c[], double fct,
	const int sign)
{
	if (plan->length == 1) return 0;
	size_t len = plan->length;
	size_t l1 = 1, nf = plan->nfct;
	cmplx* ch = RALLOC(cmplx, len);
	if (!ch) return -1;
	cmplx* p1 = c, * p2 = ch;

	for (size_t k1 = 0; k1 < nf; k1++)
	{
		size_t ip = plan->fct[k1].fct;
		size_t l2 = ip * l1;
		size_t ido = len / l2;
		if (ip == 4)
			sign > 0 ? pass4b(ido, l1, p1, p2, plan->fct[k1].tw)
			: pass4f(ido, l1, p1, p2, plan->fct[k1].tw);
		else if (ip == 2)
			sign > 0 ? pass2b(ido, l1, p1, p2, plan->fct[k1].tw)
			: pass2f(ido, l1, p1, p2, plan->fct[k1].tw);
		else if (ip == 3)
			sign > 0 ? pass3b(ido, l1, p1, p2, plan->fct[k1].tw)
			: pass3f(ido, l1, p1, p2, plan->fct[k1].tw);
		else if (ip == 5)
			sign > 0 ? pass5b(ido, l1, p1, p2, plan->fct[k1].tw)
			: pass5f(ido, l1, p1, p2, plan->fct[k1].tw);
		else if (ip == 7)  pass7(ido, l1, p1, p2, plan->fct[k1].tw, sign);
		else if (ip == 11) pass11(ido, l1, p1, p2, plan->fct[k1].tw, sign);
		else
		{
			if (passg(ido, ip, l1, p1, p2, plan->fct[k1].tw, plan->fct[k1].tws, sign))
			{
				DEALLOC(ch); return -1;
			}
			SWAP(p1, p2, cmplx*);
		}
		SWAP(p1, p2, cmplx*);
		l1 = l2;
	}
	if (p1 != c)
	{
		if (fct != 1.)
			for (size_t i = 0; i < len; ++i)
			{
				c[i].r = ch[i].r * fct;
				c[i].i = ch[i].i * fct;
			}
		else
			memcpy(c, p1, len * sizeof(cmplx));
	}
	else
		if (fct != 1.)
			for (size_t i = 0; i < len; ++i)
			{
				c[i].r *= fct;
				c[i].i *= fct;
			}
	DEALLOC(ch);
	return 0;
}

#undef PMSIGNC
#undef A_EQ_B_MUL_C
#undef A_EQ_CB_MUL_C
#undef MULPMSIGNC
#undef MULPMSIGNCEQ

#undef WA
#undef CC
#undef CH
#undef ROT90
#undef SCALEC
#undef ADDC
#undef PMC

NOINLINE WARN_UNUSED_RESULT
static int cfftp_forward(cfftp_plan plan, double c[], double fct)
{
	return pass_all(plan, (cmplx*)c, fct, -1);
}

NOINLINE WARN_UNUSED_RESULT
static int cfftp_backward(cfftp_plan plan, double c[], double fct)
{
	return pass_all(plan, (cmplx*)c, fct, 1);
}

NOINLINE WARN_UNUSED_RESULT
static int cfftp_factorize(cfftp_plan plan)
{
	size_t length = plan->length;
	size_t nfct = 0;
	while ((length % 4) == 0)
	{
		if (nfct >= NFCT) return -1; plan->fct[nfct++].fct = 4; length >>= 2;
	}
	if ((length % 2) == 0)
	{
		length >>= 1;
		// factor 2 should be at the front of the factor list
		if (nfct >= NFCT) return -1;
		plan->fct[nfct++].fct = 2;
		SWAP(plan->fct[0].fct, plan->fct[nfct - 1].fct, size_t);
	}
	size_t maxl = (size_t)(sqrt((double)length)) + 1;
	for (size_t divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
		if ((length % divisor) == 0)
		{
			while ((length % divisor) == 0)
			{
				if (nfct >= NFCT) return -1;
				plan->fct[nfct++].fct = divisor;
				length /= divisor;
			}
			maxl = (size_t)(sqrt((double)length)) + 1;
		}
	if (length > 1) plan->fct[nfct++].fct = length;
	plan->nfct = nfct;
	return 0;
}

NOINLINE static size_t cfftp_twsize(cfftp_plan plan)
{
	size_t twsize = 0, l1 = 1;
	for (size_t k = 0; k < plan->nfct; ++k)
	{
		size_t ip = plan->fct[k].fct, ido = plan->length / (l1 * ip);
		twsize += (ip - 1) * (ido - 1);
		if (ip > 11)
			twsize += ip;
		l1 *= ip;
	}
	return twsize;
}

NOINLINE WARN_UNUSED_RESULT static int cfftp_comp_twiddle(cfftp_plan plan)
{
	size_t length = plan->length;
	double* twid = RALLOC(double, 2 * length);
	if (!twid) return -1;
	sincos_2pibyn(length, twid);
	size_t l1 = 1;
	size_t memofs = 0;
	for (size_t k = 0; k < plan->nfct; ++k)
	{
		size_t ip = plan->fct[k].fct, ido = length / (l1 * ip);
		plan->fct[k].tw = plan->mem + memofs;
		memofs += (ip - 1) * (ido - 1);
		for (size_t j = 1; j < ip; ++j)
			for (size_t i = 1; i < ido; ++i)
			{
				plan->fct[k].tw[(j - 1) * (ido - 1) + i - 1].r = twid[2 * j * l1 * i];
				plan->fct[k].tw[(j - 1) * (ido - 1) + i - 1].i = twid[2 * j * l1 * i + 1];
			}
		if (ip > 11)
		{
			plan->fct[k].tws = plan->mem + memofs;
			memofs += ip;
			for (size_t j = 0; j < ip; ++j)
			{
				plan->fct[k].tws[j].r = twid[2 * j * l1 * ido];
				plan->fct[k].tws[j].i = twid[2 * j * l1 * ido + 1];
			}
		}
		l1 *= ip;
	}
	DEALLOC(twid);
	return 0;
}

static cfftp_plan make_cfftp_plan(size_t length)
{
	if (length == 0) return NULL;
	cfftp_plan plan = RALLOC(cfftp_plan_i, 1);
	if (!plan) return NULL;
	plan->length = length;
	plan->nfct = 0;
	for (size_t i = 0; i < NFCT; ++i)
		plan->fct[i] = (cfftp_fctdata){ 0,0,0 };
	plan->mem = 0;
	if (length == 1) return plan;
	if (cfftp_factorize(plan) != 0) { DEALLOC(plan); return NULL; }
	size_t tws = cfftp_twsize(plan);
	plan->mem = RALLOC(cmplx, tws);
	if (!plan->mem) { DEALLOC(plan); return NULL; }
	if (cfftp_comp_twiddle(plan) != 0)
	{
		DEALLOC(plan->mem); DEALLOC(plan); return NULL;
	}

	//if (NAIL_TEST) printf("cfftp plan mem\n");
	//if (NAIL_TEST) NailTestPrintComplexArray(plan->mem, tws);

	return plan;
}

static void destroy_cfftp_plan(cfftp_plan plan)
{
	DEALLOC(plan->mem);
	DEALLOC(plan);
}

typedef struct rfftp_fctdata
{
	size_t fct;
	double* tw, * tws;
} rfftp_fctdata;

typedef struct rfftp_plan_i
{
	size_t length, nfct;
	double* mem;
	rfftp_fctdata fct[NFCT];
} rfftp_plan_i;
typedef struct rfftp_plan_i* rfftp_plan;

#define WA(x,i) wa[(i)+(x)*(ido-1)]
#define PM(a,b,c,d) { a=c+d; b=c-d; }
/* (a+ib) = conj(c+id) * (e+if) */
#define MULPM(a,b,c,d,e,f) { a=c*e+d*f; b=c*f-d*e; }

#define CC(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radf2(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 2;

	if (NAIL_TEST) printf("start radf2\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + cc[(0) + ido * ((k)+l1 * (1))]; ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] - cc[(0) + ido * ((k)+l1 * (1))];
	}
		if ((ido & 1) == 0)
			for (size_t k = 0; k < l1; k++)
			{
				ch[(0) + ido * ((1) + cdim * (k))] = -cc[(ido - 1) + ido * ((k)+l1 * (1))];
				ch[(ido - 1) + ido * ((0) + cdim * (k))] = cc[(ido - 1) + ido * ((k)+l1 * (0))];
			}
	if (ido <= 2) return;
	for (size_t k = 0; k < l1; k++)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double tr2, ti2;
			{ tr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))]; ti2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))]; }
			{ ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k)+l1 * (0))] + tr2; ch[(ic - 1) + ido * ((1) + cdim * (k))] = cc[(i - 1) + ido * ((k)+l1 * (0))] - tr2; }
			{ ch[(i)+ido * ((0) + cdim * (k))] = ti2 + cc[(i)+ido * ((k)+l1 * (0))]; ch[(ic)+ido * ((1) + cdim * (k))] = ti2 - cc[(i)+ido * ((k)+l1 * (0))]; }
		}
}

NOINLINE static void radf3(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 3;
	static const double taur = -0.5, taui = 0.86602540378443864676;

	if (NAIL_TEST) printf("start radf3\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double cr2 = cc[(0) + ido * ((k)+l1 * (1))] + cc[(0) + ido * ((k)+l1 * (2))];
		ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + cr2;
		ch[(0) + ido * ((2) + cdim * (k))] = taui * (cc[(0) + ido * ((k)+l1 * (2))] - cc[(0) + ido * ((k)+l1 * (1))]);
		ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + taur * cr2;
	}
	if (ido == 1) return;
	for (size_t k = 0; k < l1; k++)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double di2, di3, dr2, dr3;
			{ dr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))]; di2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))]; } // d2=conj(WA0)*CC1
			{ dr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))]; di3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))]; } // d3=conj(WA1)*CC2
				double cr2 = dr2 + dr3; // c add
			double ci2 = di2 + di3;
			ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k)+l1 * (0))] + cr2; // c add
			ch[(i)+ido * ((0) + cdim * (k))] = cc[(i)+ido * ((k)+l1 * (0))] + ci2;
			double tr2 = cc[(i - 1) + ido * ((k)+l1 * (0))] + taur * cr2; // c add
			double ti2 = cc[(i)+ido * ((k)+l1 * (0))] + taur * ci2;
			double tr3 = taui * (di2 - di3);  // t3 = taui*i*(d3-d2)?
			double ti3 = taui * (dr3 - dr2);
			{ ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr3; ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr3; } // PM(i) = t2+t3
			{ ch[(i)+ido * ((2) + cdim * (k))] = ti3 + ti2; ch[(ic)+ido * ((1) + cdim * (k))] = ti3 - ti2; } // PM(ic) = conj(t2-t3)
		}
}

NOINLINE static void radf4(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 4;
	static const double hsqt2 = 0.70710678118654752440;

	if (NAIL_TEST) printf("start radf4\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double tr1, tr2;
		tr1 = cc[(0) + ido * ((k)+l1 * (3))] + cc[(0) + ido * ((k)+l1 * (1))];
		ch[(0) + ido * ((2) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (3))] - cc[(0) + ido * ((k)+l1 * (1))];
		tr2 = cc[(0) + ido * ((k)+l1 * (0))] + cc[(0) + ido * ((k)+l1 * (2))];
		ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] - cc[(0) + ido * ((k)+l1 * (2))];
		ch[(0) + ido * ((0) + cdim * (k))] = tr2 + tr1;
		ch[(ido - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
	}
	if ((ido & 1) == 0)
		for (size_t k = 0; k < l1; k++)
		{
			double ti1 = -hsqt2 * (cc[(ido - 1) + ido * ((k)+l1 * (1))] + cc[(ido - 1) + ido * ((k)+l1 * (3))]);
			double tr1 = hsqt2 * (cc[(ido - 1) + ido * ((k)+l1 * (1))] - cc[(ido - 1) + ido * ((k)+l1 * (3))]);
			ch[(ido - 1) + ido * ((0) + cdim * (k))] = cc[(ido - 1) + ido * ((k)+l1 * (0))] + tr1;
			ch[(ido - 1) + ido * ((2) + cdim * (k))] = cc[(ido - 1) + ido * ((k)+l1 * (0))] - tr1;
			ch[(0) + ido * ((3) + cdim * (k))] = ti1 + cc[(ido - 1) + ido * ((k)+l1 * (2))];
			ch[(0) + ido * ((1) + cdim * (k))] = ti1 - cc[(ido - 1) + ido * ((k)+l1 * (2))];
		}
	if (ido <= 2) return;
	for (size_t k = 0; k < l1; k++)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
			cr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))];
			ci2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))];
			cr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))];
			ci3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))];
			cr4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (3))] + wa[(i - 1) + (2) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (3))];
			ci4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (3))] - wa[(i - 1) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (3))];
			tr1 = cr4 + cr2; tr4 = cr4 - cr2;
			ti1 = ci2 + ci4; ti4 = ci2 - ci4;
			tr2 = cc[(i - 1) + ido * ((k)+l1 * (0))] + cr3;
			tr3 = cc[(i - 1) + ido * ((k)+l1 * (0))] - cr3;
			ti2 = cc[(i)+ido * ((k)+l1 * (0))] + ci3;
			ti3 = cc[(i)+ido * ((k)+l1 * (0))] - ci3;
			ch[(i - 1) + ido * ((0) + cdim * (k))] = tr2 + tr1; ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
			ch[(i)+ido * ((0) + cdim * (k))] = ti1 + ti2;
			ch[(ic)+ido * ((3) + cdim * (k))] = ti1 - ti2;
			ch[(i - 1) + ido * ((2) + cdim * (k))] = tr3 + ti4; ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr3 - ti4;
			ch[(i)+ido * ((2) + cdim * (k))] = tr4 + ti3;
			ch[(ic)+ido * ((1) + cdim * (k))] = tr4 - ti3;
		}
}

NOINLINE static void radf5(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 5;
	static const double tr11 = 0.3090169943749474241, ti11 = 0.95105651629515357212,
		tr12 = -0.8090169943749474241, ti12 = 0.58778525229247312917;

	if (NAIL_TEST) printf("start radf5\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double cr2, cr3, ci4, ci5;
		cr2 = cc[(0) + ido * ((k)+l1 * (4))] + cc[(0) + ido * ((k)+l1 * (1))];
		ci5 = cc[(0) + ido * ((k)+l1 * (4))] - cc[(0) + ido * ((k)+l1 * (1))];
		cr3 = cc[(0) + ido * ((k)+l1 * (3))] + cc[(0) + ido * ((k)+l1 * (2))];
		ci4 = cc[(0) + ido * ((k)+l1 * (3))] - cc[(0) + ido * ((k)+l1 * (2))];
		ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + cr2 + cr3;
		ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + tr11 * cr2 + tr12 * cr3;
		ch[(0) + ido * ((2) + cdim * (k))] = ti11 * ci5 + ti12 * ci4;
		ch[(ido - 1) + ido * ((3) + cdim * (k))] = cc[(0) + ido * ((k)+l1 * (0))] + tr12 * cr2 + tr11 * cr3;
		ch[(0) + ido * ((4) + cdim * (k))] = ti12 * ci5 - ti11 * ci4;
	}
	if (ido == 1) return;
	for (size_t k = 0; k < l1; ++k)
		for (size_t i = 2; i < ido; i += 2)
		{
			double ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3,
				dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
			size_t ic = ido - i;
			dr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))];
			di2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (1))];
			dr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))];
			di3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (2))];
			dr4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (3))] + wa[(i - 1) + (2) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (3))];
			di4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (3))] - wa[(i - 1) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (3))];
			dr5 = wa[(i - 2) + (3) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (4))] + wa[(i - 1) + (3) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (4))];
			di5 = wa[(i - 2) + (3) * (ido - 1)] * cc[(i)+ido * ((k)+l1 * (4))] - wa[(i - 1) + (3) * (ido - 1)] * cc[(i - 1) + ido * ((k)+l1 * (4))];
			cr2 = dr5 + dr2; ci5 = dr5 - dr2;
			ci2 = di2 + di5; cr5 = di2 - di5;
			cr3 = dr4 + dr3; ci4 = dr4 - dr3;
			ci3 = di3 + di4; cr4 = di3 - di4;
			ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k)+l1 * (0))] + cr2 + cr3;
			ch[(i)+ido * ((0) + cdim * (k))] = cc[(i)+ido * ((k)+l1 * (0))] + ci2 + ci3;
			tr2 = cc[(i - 1) + ido * ((k)+l1 * (0))] + tr11 * cr2 + tr12 * cr3;
			ti2 = cc[(i)+ido * ((k)+l1 * (0))] + tr11 * ci2 + tr12 * ci3;
			tr3 = cc[(i - 1) + ido * ((k)+l1 * (0))] + tr12 * cr2 + tr11 * cr3;
			ti3 = cc[(i)+ido * ((k)+l1 * (0))] + tr12 * ci2 + tr11 * ci3;
			tr5 = cr5 * ti11 + cr4 * ti12;
			tr4 = cr5 * ti12 - cr4 * ti11;
			ti5 = ci5 * ti11 + ci4 * ti12;
			ti4 = ci5 * ti12 - ci4 * ti11;
			ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr5;
			ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr5;
			ch[(i)+ido * ((2) + cdim * (k))] = ti5 + ti2;
			ch[(ic)+ido * ((1) + cdim * (k))] = ti5 - ti2;
			ch[(i - 1) + ido * ((4) + cdim * (k))] = tr3 + tr4;
			ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr3 - tr4;
			ch[(i)+ido * ((4) + cdim * (k))] = ti4 + ti3;
			ch[(ic)+ido * ((3) + cdim * (k))] = ti4 - ti3;
		}
}

#undef CC
#undef CH
#define C1(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define C2(a,b) cc[(a)+idl1*(b)]
#define CH2(a,b) ch[(a)+idl1*(b)]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
NOINLINE static void radfg(size_t ido, size_t ip, size_t l1,
	double* cc, double* ch, const double* wa,
	const double* csarr)
{
	const size_t cdim = ip;
	size_t ipph = (ip + 1) / 2;
	size_t idl1 = ido * l1;

	if (NAIL_TEST) printf("start radfg\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	if (ido > 1)
	{
		for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)              // 114
		{
			size_t is = (j - 1) * (ido - 1),
				is2 = (jc - 1) * (ido - 1);
			for (size_t k = 0; k < l1; ++k)                            // 113
			{
				size_t idij = is;
				size_t idij2 = is2;
				for (size_t i = 1; i <= ido - 2; i += 2)                      // 112
				{
					double t1 = cc[(i)+ido * ((k)+l1 * (j))], t2 = cc[(i + 1) + ido * ((k)+l1 * (j))],
						t3 = cc[(i)+ido * ((k)+l1 * (jc))], t4 = cc[(i + 1) + ido * ((k)+l1 * (jc))];
					double x1 = wa[idij] * t1 + wa[idij + 1] * t2,
						x2 = wa[idij] * t2 - wa[idij + 1] * t1,
						x3 = wa[idij2] * t3 + wa[idij2 + 1] * t4,
						x4 = wa[idij2] * t4 - wa[idij2 + 1] * t3;
					cc[(i)+ido * ((k)+l1 * (j))] = x1 + x3;
					cc[(i)+ido * ((k)+l1 * (jc))] = x2 - x4;
					cc[(i + 1) + ido * ((k)+l1 * (j))] = x2 + x4;
					cc[(i + 1) + ido * ((k)+l1 * (jc))] = x3 - x1;
					idij += 2;
					idij2 += 2;
				}
			}
		}
	}

	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 123
		for (size_t k = 0; k < l1; ++k)                              // 122
		{
			double t1 = cc[(0) + ido * ((k)+l1 * (j))], t2 = cc[(0) + ido * ((k)+l1 * (jc))];
			cc[(0) + ido * ((k)+l1 * (j))] = t1 + t2;
			cc[(0) + ido * ((k)+l1 * (jc))] = t2 - t1;
		}

	//everything in C
	//memset(ch,0,ip*l1*ido*sizeof(double));

	for (size_t l = 1, lc = ip - 1; l < ipph; ++l, --lc)                 // 127
	{
		for (size_t ik = 0; ik < idl1; ++ik)                         // 124
		{
			ch[(ik)+idl1 * (l)] = cc[(ik)+idl1 * (0)] + csarr[2 * l] * cc[(ik)+idl1 * (1)] + csarr[4 * l] * cc[(ik)+idl1 * (2)];
			ch[(ik)+idl1 * (lc)] = csarr[2 * l + 1] * cc[(ik)+idl1 * (ip - 1)] + csarr[4 * l + 1] * cc[(ik)+idl1 * (ip - 2)];
		}

		size_t iang = 2 * l;
		size_t j = 3, jc = ip - 3;
		for (; j < ipph - 3; j += 4, jc -= 4)              // 126
		{
			iang += l; if (iang >= ip) iang -= ip;
			double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
			iang += l; if (iang >= ip) iang -= ip;
			double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
			iang += l; if (iang >= ip) iang -= ip;
			double ar3 = csarr[2 * iang], ai3 = csarr[2 * iang + 1];
			iang += l; if (iang >= ip) iang -= ip;
			double ar4 = csarr[2 * iang], ai4 = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)                       // 125
			{
				ch[(ik)+idl1 * (l)] += ar1 * cc[(ik)+idl1 * (j)] + ar2 * cc[(ik)+idl1 * (j + 1)]
					+ ar3 * cc[(ik)+idl1 * (j + 2)] + ar4 * cc[(ik)+idl1 * (j + 3)];
				ch[(ik)+idl1 * (lc)] += ai1 * cc[(ik)+idl1 * (jc)] + ai2 * cc[(ik)+idl1 * (jc - 1)]
					+ ai3 * cc[(ik)+idl1 * (jc - 2)] + ai4 * cc[(ik)+idl1 * (jc - 3)];
			}
		}
		for (; j < ipph - 1; j += 2, jc -= 2)              // 126
		{
			iang += l; if (iang >= ip) iang -= ip;
			double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
			iang += l; if (iang >= ip) iang -= ip;
			double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)                       // 125
			{
				ch[(ik)+idl1 * (l)] += ar1 * cc[(ik)+idl1 * (j)] + ar2 * cc[(ik)+idl1 * (j + 1)];
				ch[(ik)+idl1 * (lc)] += ai1 * cc[(ik)+idl1 * (jc)] + ai2 * cc[(ik)+idl1 * (jc - 1)];
			}
		}
		for (; j < ipph; ++j, --jc)              // 126
		{
			iang += l; if (iang >= ip) iang -= ip;
			double ar = csarr[2 * iang], ai = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)                       // 125
			{
				ch[(ik)+idl1 * (l)] += ar * cc[(ik)+idl1 * (j)];
				ch[(ik)+idl1 * (lc)] += ai * cc[(ik)+idl1 * (jc)];
			}
		}
	}
	for (size_t ik = 0; ik < idl1; ++ik)                         // 101
		ch[(ik)+idl1 * (0)] = cc[(ik)+idl1 * (0)];
	for (size_t j = 1; j < ipph; ++j)                              // 129
		for (size_t ik = 0; ik < idl1; ++ik)                         // 128
			ch[(ik)+idl1 * (0)] += cc[(ik)+idl1 * (j)];

	// everything in CH at this point!
	//memset(cc,0,ip*l1*ido*sizeof(double));

	for (size_t k = 0; k < l1; ++k)                                // 131
		for (size_t i = 0; i < ido; ++i)                             // 130
			cc[(i)+ido * ((0) + cdim * (k))] = ch[(i)+ido * ((k)+l1 * (0))];

	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 137
	{
		size_t j2 = 2 * j - 1;
		for (size_t k = 0; k < l1; ++k)                              // 136
		{
			cc[(ido - 1) + ido * ((j2)+cdim * (k))] = ch[(0) + ido * ((k)+l1 * (j))];
			cc[(0) + ido * ((j2 + 1) + cdim * (k))] = ch[(0) + ido * ((k)+l1 * (jc))];
		}
	}

	if (ido == 1) return;

	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 140
	{
		size_t j2 = 2 * j - 1;
		for (size_t k = 0; k < l1; ++k)                               // 139
			for (size_t i = 1, ic = ido - i - 2; i <= ido - 2; i += 2, ic -= 2)      // 138
			{
				cc[(i)+ido * ((j2 + 1) + cdim * (k))] = ch[(i)+ido * ((k)+l1 * (j))] + ch[(i)+ido * ((k)+l1 * (jc))];
				cc[(ic)+ido * ((j2)+cdim * (k))] = ch[(i)+ido * ((k)+l1 * (j))] - ch[(i)+ido * ((k)+l1 * (jc))];
				cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] = ch[(i + 1) + ido * ((k)+l1 * (j))] + ch[(i + 1) + ido * ((k)+l1 * (jc))];
				cc[(ic + 1) + ido * ((j2)+cdim * (k))] = ch[(i + 1) + ido * ((k)+l1 * (jc))] - ch[(i + 1) + ido * ((k)+l1 * (j))];
			}
	}
}
#undef C1
#undef C2
#undef CH2

#undef CH
#undef CC
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]

NOINLINE static void radb2(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 2;

	if (NAIL_TEST) printf("start radb2\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		ch[(0) + ido * ((k)+l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((1) + cdim * (k))]; ch[(0) + ido * ((k)+l1 * (1))] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((1) + cdim * (k))];
	}
		if ((ido & 1) == 0)
			for (size_t k = 0; k < l1; k++)
			{
				ch[(ido - 1) + ido * ((k)+l1 * (0))] = 2. * cc[(ido - 1) + ido * ((0) + cdim * (k))];
				ch[(ido - 1) + ido * ((k)+l1 * (1))] = -2. * cc[(0) + ido * ((1) + cdim * (k))];
			}
	if (ido <= 2) return;
	for (size_t k = 0; k < l1; ++k)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double ti2, tr2;
			{ ch[(i - 1) + ido * ((k)+l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; tr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))]; }
			{ ti2 = cc[(i)+ido * ((0) + cdim * (k))] + cc[(ic)+ido * ((1) + cdim * (k))]; ch[(i)+ido * ((k)+l1 * (0))] = cc[(i)+ido * ((0) + cdim * (k))] - cc[(ic)+ido * ((1) + cdim * (k))]; }
			{ ch[(i)+ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * ti2 + wa[(i - 1) + (0) * (ido - 1)] * tr2; ch[(i - 1) + ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * tr2 - wa[(i - 1) + (0) * (ido - 1)] * ti2; }
		}
}

NOINLINE static void radb3(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 3;
	static const double taur = -0.5, taui = 0.86602540378443864676;

	if (NAIL_TEST) printf("start radb3\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double tr2 = 2. * cc[(ido - 1) + ido * ((1) + cdim * (k))];
		double cr2 = cc[(0) + ido * ((0) + cdim * (k))] + taur * tr2;
		ch[(0) + ido * ((k)+l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + tr2;
		double ci3 = 2. * taui * cc[(0) + ido * ((2) + cdim * (k))];
		{ ch[(0) + ido * ((k)+l1 * (2))] = cr2 + ci3; ch[(0) + ido * ((k)+l1 * (1))] = cr2 - ci3; };
	}
	if (ido == 1) return;
	for (size_t k = 0; k < l1; k++)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double tr2 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; // t2=CC(I) + conj(CC(ic))
			double ti2 = cc[(i)+ido * ((2) + cdim * (k))] - cc[(ic)+ido * ((1) + cdim * (k))];
			double cr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + taur * tr2;     // c2=CC +taur*t2
			double ci2 = cc[(i)+ido * ((0) + cdim * (k))] + taur * ti2;
			ch[(i - 1) + ido * ((k)+l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2;         // CH=CC+t2
			ch[(i)+ido * ((k)+l1 * (0))] = cc[(i)+ido * ((0) + cdim * (k))] + ti2;
			double cr3 = taui * (cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))]);// c3=taui*(CC(i)-conj(CC(ic)))
			double ci3 = taui * (cc[(i)+ido * ((2) + cdim * (k))] + cc[(ic)+ido * ((1) + cdim * (k))]);
			double di2, di3, dr2, dr3;
			{ dr3 = cr2 + ci3; dr2 = cr2 - ci3; } // d2= (cr2-ci3, ci2+cr3) = c2+i*c3
			{ di2 = ci2 + cr3; di3 = ci2 - cr3; } // d3= (cr2+ci3, ci2-cr3) = c2-i*c3
			{ ch[(i)+ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * di2 + wa[(i - 1) + (0) * (ido - 1)] * dr2; ch[(i - 1) + ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * dr2 - wa[(i - 1) + (0) * (ido - 1)] * di2; } // ch = WA*d2
			{ ch[(i)+ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * di3 + wa[(i - 1) + (1) * (ido - 1)] * dr3; ch[(i - 1) + ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * dr3 - wa[(i - 1) + (1) * (ido - 1)] * di3; }
		}
}

NOINLINE static void radb4(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 4;
	static const double sqrt2 = 1.41421356237309504880;

	if (NAIL_TEST) printf("start radb4\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double tr1, tr2;
		{ tr2 = cc[(0) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((3) + cdim * (k))]; tr1 = cc[(0) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((3) + cdim * (k))]; }
		double tr3 = 2. * cc[(ido - 1) + ido * ((1) + cdim * (k))];
		double tr4 = 2. * cc[(0) + ido * ((2) + cdim * (k))];
		{ ch[(0) + ido * ((k)+l1 * (0))] = tr2 + tr3; ch[(0) + ido * ((k)+l1 * (2))] = tr2 - tr3; }
		{ ch[(0) + ido * ((k)+l1 * (3))] = tr1 + tr4; ch[(0) + ido * ((k)+l1 * (1))] = tr1 - tr4; }
	}
	if ((ido & 1) == 0)
		for (size_t k = 0; k < l1; k++)
		{
			double tr1, tr2, ti1, ti2;
			{ ti1 = cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((1) + cdim * (k))]; ti2 = cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))]; }
			{ tr2 = cc[(ido - 1) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((2) + cdim * (k))]; tr1 = cc[(ido - 1) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((2) + cdim * (k))]; }
			ch[(ido - 1) + ido * ((k)+l1 * (0))] = tr2 + tr2;
			ch[(ido - 1) + ido * ((k)+l1 * (1))] = sqrt2 * (tr1 - ti1);
			ch[(ido - 1) + ido * ((k)+l1 * (2))] = ti2 + ti2;
			ch[(ido - 1) + ido * ((k)+l1 * (3))] = -sqrt2 * (tr1 + ti1);
		}
	if (ido <= 2) return;
	for (size_t k = 0; k < l1; ++k)
		for (size_t i = 2; i < ido; i += 2)
		{
			double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
			size_t ic = ido - i;
			{ tr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + cc[(ic - 1) + ido * ((3) + cdim * (k))]; tr1 = cc[(i - 1) + ido * ((0) + cdim * (k))] - cc[(ic - 1) + ido * ((3) + cdim * (k))]; }
			{ ti1 = cc[(i)+ido * ((0) + cdim * (k))] + cc[(ic)+ido * ((3) + cdim * (k))]; ti2 = cc[(i)+ido * ((0) + cdim * (k))] - cc[(ic)+ido * ((3) + cdim * (k))]; }
			{ tr4 = cc[(i)+ido * ((2) + cdim * (k))] + cc[(ic)+ido * ((1) + cdim * (k))]; ti3 = cc[(i)+ido * ((2) + cdim * (k))] - cc[(ic)+ido * ((1) + cdim * (k))]; }
			{ tr3 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; ti4 = cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))]; }
			{ ch[(i - 1) + ido * ((k)+l1 * (0))] = tr2 + tr3; cr3 = tr2 - tr3; }
			{ ch[(i)+ido * ((k)+l1 * (0))] = ti2 + ti3; ci3 = ti2 - ti3; }
			{ cr4 = tr1 + tr4; cr2 = tr1 - tr4; }
			{ ci2 = ti1 + ti4; ci4 = ti1 - ti4; }
			{ ch[(i)+ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * ci2 + wa[(i - 1) + (0) * (ido - 1)] * cr2; ch[(i - 1) + ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * cr2 - wa[(i - 1) + (0) * (ido - 1)] * ci2; }
			{ ch[(i)+ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * ci3 + wa[(i - 1) + (1) * (ido - 1)] * cr3; ch[(i - 1) + ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * cr3 - wa[(i - 1) + (1) * (ido - 1)] * ci3; }
			{ ch[(i)+ido * ((k)+l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * ci4 + wa[(i - 1) + (2) * (ido - 1)] * cr4; ch[(i - 1) + ido * ((k)+l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * cr4 - wa[(i - 1) + (2) * (ido - 1)] * ci4; }
		}
}

NOINLINE static void radb5(size_t ido, size_t l1, const double* cc,
	double* ch, const double* wa)
{
	const size_t cdim = 5;
	static const double tr11 = 0.3090169943749474241, ti11 = 0.95105651629515357212,
		tr12 = -0.8090169943749474241, ti12 = 0.58778525229247312917;

	if (NAIL_TEST) printf("start radb5\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; k++)
	{
		double ti5 = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
		double ti4 = cc[(0) + ido * ((4) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
		double tr2 = cc[(ido - 1) + ido * ((1) + cdim * (k))] + cc[(ido - 1) + ido * ((1) + cdim * (k))];
		double tr3 = cc[(ido - 1) + ido * ((3) + cdim * (k))] + cc[(ido - 1) + ido * ((3) + cdim * (k))];
		ch[(0) + ido * ((k)+l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + tr2 + tr3;
		double cr2 = cc[(0) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
		double cr3 = cc[(0) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
		double ci4, ci5;
		{ ci5 = ti5 * ti11 + ti4 * ti12; ci4 = ti5 * ti12 - ti4 * ti11; }
		{ ch[(0) + ido * ((k)+l1 * (4))] = cr2 + ci5; ch[(0) + ido * ((k)+l1 * (1))] = cr2 - ci5; }
		{ ch[(0) + ido * ((k)+l1 * (3))] = cr3 + ci4; ch[(0) + ido * ((k)+l1 * (2))] = cr3 - ci4; }
	}
	if (ido == 1) return;
	for (size_t k = 0; k < l1; ++k)
		for (size_t i = 2; i < ido; i += 2)
		{
			size_t ic = ido - i;
			double tr2, tr3, tr4, tr5, ti2, ti3, ti4, ti5;
			{ tr2 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; tr5 = cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))]; }
			{ ti5 = cc[(i)+ido * ((2) + cdim * (k))] + cc[(ic)+ido * ((1) + cdim * (k))]; ti2 = cc[(i)+ido * ((2) + cdim * (k))] - cc[(ic)+ido * ((1) + cdim * (k))]; }
			{ tr3 = cc[(i - 1) + ido * ((4) + cdim * (k))] + cc[(ic - 1) + ido * ((3) + cdim * (k))]; tr4 = cc[(i - 1) + ido * ((4) + cdim * (k))] - cc[(ic - 1) + ido * ((3) + cdim * (k))]; }
			{ ti4 = cc[(i)+ido * ((4) + cdim * (k))] + cc[(ic)+ido * ((3) + cdim * (k))]; ti3 = cc[(i)+ido * ((4) + cdim * (k))] - cc[(ic)+ido * ((3) + cdim * (k))]; }
			ch[(i - 1) + ido * ((k)+l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2 + tr3;
			ch[(i)+ido * ((k)+l1 * (0))] = cc[(i)+ido * ((0) + cdim * (k))] + ti2 + ti3;
			double cr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
			double ci2 = cc[(i)+ido * ((0) + cdim * (k))] + tr11 * ti2 + tr12 * ti3;
			double cr3 = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
			double ci3 = cc[(i)+ido * ((0) + cdim * (k))] + tr12 * ti2 + tr11 * ti3;
			double ci4, ci5, cr5, cr4;
			{ cr5 = tr5 * ti11 + tr4 * ti12; cr4 = tr5 * ti12 - tr4 * ti11; }
			{ ci5 = ti5 * ti11 + ti4 * ti12; ci4 = ti5 * ti12 - ti4 * ti11; }
				double dr2, dr3, dr4, dr5, di2, di3, di4, di5;
				dr4 = cr3 + ci4; dr3 = cr3 - ci4;
				di3 = ci3 + cr4; di4 = ci3 - cr4;
				dr5 = cr2 + ci5; dr2 = cr2 - ci5;
				di2 = ci2 + cr5; di5 = ci2 - cr5;
				ch[(i)+ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * di2 + wa[(i - 1) + (0) * (ido - 1)] * dr2;
				ch[(i - 1) + ido * ((k)+l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * dr2 - wa[(i - 1) + (0) * (ido - 1)] * di2;
				ch[(i)+ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * di3 + wa[(i - 1) + (1) * (ido - 1)] * dr3;
				ch[(i - 1) + ido * ((k)+l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * dr3 - wa[(i - 1) + (1) * (ido - 1)] * di3;
				ch[(i)+ido * ((k)+l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * di4 + wa[(i - 1) + (2) * (ido - 1)] * dr4;
				ch[(i - 1) + ido * ((k)+l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * dr4 - wa[(i - 1) + (2) * (ido - 1)] * di4;
				ch[(i)+ido * ((k)+l1 * (4))] = wa[(i - 2) + (3) * (ido - 1)] * di5 + wa[(i - 1) + (3) * (ido - 1)] * dr5;
				ch[(i - 1) + ido * ((k)+l1 * (4))] = wa[(i - 2) + (3) * (ido - 1)] * dr5 - wa[(i - 1) + (3) * (ido - 1)] * di5;
		}
}

#undef CC
#undef CH
#define CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
#define CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
#define C1(a,b,c) cc[(a)+ido*((b)+l1*(c))]
#define C2(a,b) cc[(a)+idl1*(b)]
#define CH2(a,b) ch[(a)+idl1*(b)]

NOINLINE static void radbg(size_t ido, size_t ip, size_t l1,
	double* cc, double* ch, const double* wa,
	const double* csarr)
{
	const size_t cdim = ip;
	size_t ipph = (ip + 1) / 2;
	size_t idl1 = ido * l1;

	if (NAIL_TEST) printf("start radbg\n");
	if (NAIL_TEST) NailTestPrintDoubleArray(cc, cdim * l1);
	for (size_t k = 0; k < l1; ++k)        // 102
		for (size_t i = 0; i < ido; ++i)     // 101
			ch[(i)+ido * ((k)+l1 * (0))] = cc[(i)+ido * ((0) + cdim * (k))];
	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 108
	{
		size_t j2 = 2 * j - 1;
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (j))] = 2 * cc[(ido - 1) + ido * ((j2)+cdim * (k))];
			ch[(0) + ido * ((k)+l1 * (jc))] = 2 * cc[(0) + ido * ((j2 + 1) + cdim * (k))];
		}
	}

	if (ido != 1)
	{
		for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 111
		{
			size_t j2 = 2 * j - 1;
			for (size_t k = 0; k < l1; ++k)
				for (size_t i = 1, ic = ido - i - 2; i <= ido - 2; i += 2, ic -= 2)      // 109
				{
					ch[(i)+ido * ((k)+l1 * (j))] = cc[(i)+ido * ((j2 + 1) + cdim * (k))] + cc[(ic)+ido * ((j2)+cdim * (k))];
					ch[(i)+ido * ((k)+l1 * (jc))] = cc[(i)+ido * ((j2 + 1) + cdim * (k))] - cc[(ic)+ido * ((j2)+cdim * (k))];
					ch[(i + 1) + ido * ((k)+l1 * (j))] = cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] - cc[(ic + 1) + ido * ((j2)+cdim * (k))];
					ch[(i + 1) + ido * ((k)+l1 * (jc))] = cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] + cc[(ic + 1) + ido * ((j2)+cdim * (k))];
				}
		}
	}
	for (size_t l = 1, lc = ip - 1; l < ipph; ++l, --lc)
	{
		for (size_t ik = 0; ik < idl1; ++ik)
		{
			cc[(ik)+idl1 * (l)] = ch[(ik)+idl1 * (0)] + csarr[2 * l] * ch[(ik)+idl1 * (1)] + csarr[4 * l] * ch[(ik)+idl1 * (2)];
			cc[(ik)+idl1 * (lc)] = csarr[2 * l + 1] * ch[(ik)+idl1 * (ip - 1)] + csarr[4 * l + 1] * ch[(ik)+idl1 * (ip - 2)];
		}
		size_t iang = 2 * l;
		size_t j = 3, jc = ip - 3;
		for (; j < ipph - 3; j += 4, jc -= 4)
		{
			iang += l; if (iang > ip) iang -= ip;
			double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
			iang += l; if (iang > ip) iang -= ip;
			double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
			iang += l; if (iang > ip) iang -= ip;
			double ar3 = csarr[2 * iang], ai3 = csarr[2 * iang + 1];
			iang += l; if (iang > ip) iang -= ip;
			double ar4 = csarr[2 * iang], ai4 = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cc[(ik)+idl1 * (l)] += ar1 * ch[(ik)+idl1 * (j)] + ar2 * ch[(ik)+idl1 * (j + 1)]
					+ ar3 * ch[(ik)+idl1 * (j + 2)] + ar4 * ch[(ik)+idl1 * (j + 3)];
				cc[(ik)+idl1 * (lc)] += ai1 * ch[(ik)+idl1 * (jc)] + ai2 * ch[(ik)+idl1 * (jc - 1)]
					+ ai3 * ch[(ik)+idl1 * (jc - 2)] + ai4 * ch[(ik)+idl1 * (jc - 3)];
			}
		}
		for (; j < ipph - 1; j += 2, jc -= 2)
		{
			iang += l; if (iang > ip) iang -= ip;
			double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
			iang += l; if (iang > ip) iang -= ip;
			double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cc[(ik)+idl1 * (l)] += ar1 * ch[(ik)+idl1 * (j)] + ar2 * ch[(ik)+idl1 * (j + 1)];
				cc[(ik)+idl1 * (lc)] += ai1 * ch[(ik)+idl1 * (jc)] + ai2 * ch[(ik)+idl1 * (jc - 1)];
			}
		}
		for (; j < ipph; ++j, --jc)
		{
			iang += l; if (iang > ip) iang -= ip;
			double war = csarr[2 * iang], wai = csarr[2 * iang + 1];
			for (size_t ik = 0; ik < idl1; ++ik)
			{
				cc[(ik)+idl1 * (l)] += war * ch[(ik)+idl1 * (j)];
				cc[(ik)+idl1 * (lc)] += wai * ch[(ik)+idl1 * (jc)];
			}
		}
	}
	for (size_t j = 1; j < ipph; ++j)
		for (size_t ik = 0; ik < idl1; ++ik)
			ch[(ik)+idl1 * (0)] += ch[(ik)+idl1 * (j)];
	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 124
		for (size_t k = 0; k < l1; ++k)
		{
			ch[(0) + ido * ((k)+l1 * (j))] = cc[(0) + ido * ((k)+l1 * (j))] - cc[(0) + ido * ((k)+l1 * (jc))];
			ch[(0) + ido * ((k)+l1 * (jc))] = cc[(0) + ido * ((k)+l1 * (j))] + cc[(0) + ido * ((k)+l1 * (jc))];
		}

	if (ido == 1) return;

	for (size_t j = 1, jc = ip - 1; j < ipph; ++j, --jc)  // 127
		for (size_t k = 0; k < l1; ++k)
			for (size_t i = 1; i <= ido - 2; i += 2)
			{
				ch[(i)+ido * ((k)+l1 * (j))] = cc[(i)+ido * ((k)+l1 * (j))] - cc[(i + 1) + ido * ((k)+l1 * (jc))];
				ch[(i)+ido * ((k)+l1 * (jc))] = cc[(i)+ido * ((k)+l1 * (j))] + cc[(i + 1) + ido * ((k)+l1 * (jc))];
				ch[(i + 1) + ido * ((k)+l1 * (j))] = cc[(i + 1) + ido * ((k)+l1 * (j))] + cc[(i)+ido * ((k)+l1 * (jc))];
				ch[(i + 1) + ido * ((k)+l1 * (jc))] = cc[(i + 1) + ido * ((k)+l1 * (j))] - cc[(i)+ido * ((k)+l1 * (jc))];
			}

	// All in CH

	for (size_t j = 1; j < ip; ++j)
	{
		size_t is = (j - 1) * (ido - 1);
		for (size_t k = 0; k < l1; ++k)
		{
			size_t idij = is;
			for (size_t i = 1; i <= ido - 2; i += 2)
			{
				double t1 = ch[(i)+ido * ((k)+l1 * (j))], t2 = ch[(i + 1) + ido * ((k)+l1 * (j))];
				ch[(i)+ido * ((k)+l1 * (j))] = wa[idij] * t1 - wa[idij + 1] * t2;
				ch[(i + 1) + ido * ((k)+l1 * (j))] = wa[idij] * t2 + wa[idij + 1] * t1;
				idij += 2;
			}
		}
	}
}
#undef C1
#undef C2
#undef CH2

#undef CC
#undef CH
#undef PM
#undef MULPM
#undef WA

static void copy_and_norm(double* c, double* p1, size_t n, double fct)
{
	if (p1 != c)
	{
		if (fct != 1.)
			for (size_t i = 0; i < n; ++i)
				c[i] = fct * p1[i];
		else
			memcpy(c, p1, n * sizeof(double));
	}
	else
		if (fct != 1.)
			for (size_t i = 0; i < n; ++i)
				c[i] *= fct;
}

WARN_UNUSED_RESULT
static int rfftp_forward(rfftp_plan plan, double c[], double fct)
{
	if (plan->length == 1) return 0;
	size_t n = plan->length;
	size_t l1 = n, nf = plan->nfct;
	double* ch = RALLOC(double, n);
	if (!ch) return -1;
	double* p1 = c, * p2 = ch;

	for (size_t k1 = 0; k1 < nf; ++k1)
	{
		size_t k = nf - k1 - 1;
		size_t ip = plan->fct[k].fct;
		size_t ido = n / l1;
		l1 /= ip;
		if (ip == 4)
			radf4(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 2)
			radf2(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 3)
			radf3(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 5)
			radf5(ido, l1, p1, p2, plan->fct[k].tw);
		else
		{
			radfg(ido, ip, l1, p1, p2, plan->fct[k].tw, plan->fct[k].tws);
			SWAP(p1, p2, double*);
		}
		SWAP(p1, p2, double*);
	}
	copy_and_norm(c, p1, n, fct);
	DEALLOC(ch);
	return 0;
}

WARN_UNUSED_RESULT
static int rfftp_backward(rfftp_plan plan, double c[], double fct)
{
	if (plan->length == 1) return 0;
	size_t n = plan->length;
	size_t l1 = 1, nf = plan->nfct;
	double* ch = RALLOC(double, n);
	if (!ch) return -1;
	double* p1 = c, * p2 = ch;

	for (size_t k = 0; k < nf; k++)
	{
		size_t ip = plan->fct[k].fct,
			ido = n / (ip * l1);
		if (ip == 4)
			radb4(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 2)
			radb2(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 3)
			radb3(ido, l1, p1, p2, plan->fct[k].tw);
		else if (ip == 5)
			radb5(ido, l1, p1, p2, plan->fct[k].tw);
		else
			radbg(ido, ip, l1, p1, p2, plan->fct[k].tw, plan->fct[k].tws);
		SWAP(p1, p2, double*);
		l1 *= ip;
	}
	copy_and_norm(c, p1, n, fct);
	DEALLOC(ch);
	return 0;
}

WARN_UNUSED_RESULT
static int rfftp_factorize(rfftp_plan plan)
{
	size_t length = plan->length;
	size_t nfct = 0;
	while ((length % 4) == 0)
	{
		if (nfct >= NFCT) return -1; plan->fct[nfct++].fct = 4; length >>= 2;
	}
	if ((length % 2) == 0)
	{
		length >>= 1;
		// factor 2 should be at the front of the factor list
		if (nfct >= NFCT) return -1;
		plan->fct[nfct++].fct = 2;
		SWAP(plan->fct[0].fct, plan->fct[nfct - 1].fct, size_t);
	}
	size_t maxl = (size_t)(sqrt((double)length)) + 1;
	for (size_t divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
		if ((length % divisor) == 0)
		{
			while ((length % divisor) == 0)
			{
				if (nfct >= NFCT) return -1;
				plan->fct[nfct++].fct = divisor;
				length /= divisor;
			}
			maxl = (size_t)(sqrt((double)length)) + 1;
		}
	if (length > 1) plan->fct[nfct++].fct = length;
	plan->nfct = nfct;
	return 0;
}

static size_t rfftp_twsize(rfftp_plan plan)
{
	size_t twsize = 0, l1 = 1;
	for (size_t k = 0; k < plan->nfct; ++k)
	{
		size_t ip = plan->fct[k].fct, ido = plan->length / (l1 * ip);
		twsize += (ip - 1) * (ido - 1);
		if (ip > 5) twsize += 2 * ip;
		l1 *= ip;
	}
	return twsize;
	return 0;
}

WARN_UNUSED_RESULT NOINLINE static int rfftp_comp_twiddle(rfftp_plan plan)
{
	size_t length = plan->length;
	double* twid = RALLOC(double, 2 * length);
	if (!twid) return -1;
	sincos_2pibyn_half(length, twid);
	size_t l1 = 1;
	double* ptr = plan->mem;
	for (size_t k = 0; k < plan->nfct; ++k)
	{
		size_t ip = plan->fct[k].fct, ido = length / (l1 * ip);
		if (k < plan->nfct - 1) // last factor doesn't need twiddles
		{
			plan->fct[k].tw = ptr; ptr += (ip - 1) * (ido - 1);
			for (size_t j = 1; j < ip; ++j)
				for (size_t i = 1; i <= (ido - 1) / 2; ++i)
				{
					plan->fct[k].tw[(j - 1) * (ido - 1) + 2 * i - 2] = twid[2 * j * l1 * i];
					plan->fct[k].tw[(j - 1) * (ido - 1) + 2 * i - 1] = twid[2 * j * l1 * i + 1];
				}
		}
		if (ip > 5) // special factors required by *g functions
		{
			plan->fct[k].tws = ptr; ptr += 2 * ip;
			plan->fct[k].tws[0] = 1.;
			plan->fct[k].tws[1] = 0.;
			for (size_t i = 1; i <= (ip >> 1); ++i)
			{
				plan->fct[k].tws[2 * i] = twid[2 * i * (length / ip)];
				plan->fct[k].tws[2 * i + 1] = twid[2 * i * (length / ip) + 1];
				plan->fct[k].tws[2 * (ip - i)] = twid[2 * i * (length / ip)];
				plan->fct[k].tws[2 * (ip - i) + 1] = -twid[2 * i * (length / ip) + 1];
			}
		}
		l1 *= ip;
	}
	DEALLOC(twid);
	return 0;
}

NOINLINE static rfftp_plan make_rfftp_plan(size_t length)
{
	if (length == 0) return NULL;
	rfftp_plan plan = RALLOC(rfftp_plan_i, 1);
	if (!plan) return NULL;
	plan->length = length;
	plan->nfct = 0;
	plan->mem = NULL;
	for (size_t i = 0; i < NFCT; ++i)
		plan->fct[i] = (rfftp_fctdata){ 0,0,0 };
	if (length == 1) return plan;
	if (rfftp_factorize(plan) != 0) { DEALLOC(plan); return NULL; }
	size_t tws = rfftp_twsize(plan);
	plan->mem = RALLOC(double, tws);
	if (!plan->mem) { DEALLOC(plan); return NULL; }
	if (rfftp_comp_twiddle(plan) != 0)
	{
		DEALLOC(plan->mem); DEALLOC(plan); return NULL;
	}

	//if (NAIL_TEST) printf("rfftp plan mem\n");
	//if (NAIL_TEST) NailTestPrintDoubleArray(plan->mem, tws);

	return plan;
}

NOINLINE static void destroy_rfftp_plan(rfftp_plan plan)
{
	DEALLOC(plan->mem);
	DEALLOC(plan);
}

typedef struct fftblue_plan_i
{
	size_t n, n2;
	cfftp_plan plan;
	double* mem;
	double* bk, * bkf;
} fftblue_plan_i;
typedef struct fftblue_plan_i* fftblue_plan;

NOINLINE static fftblue_plan make_fftblue_plan(size_t length)
{
	fftblue_plan plan = RALLOC(fftblue_plan_i, 1);
	if (!plan) return NULL;
	plan->n = length;
	plan->n2 = good_size(plan->n * 2 - 1);
	plan->mem = RALLOC(double, 2 * plan->n + 2 * plan->n2);
	if (!plan->mem) { DEALLOC(plan); return NULL; }
	plan->bk = plan->mem;
	plan->bkf = plan->bk + 2 * plan->n;

	/* initialize b_k */
	double* tmp = RALLOC(double, 4 * plan->n);
	if (!tmp) { DEALLOC(plan->mem); DEALLOC(plan); return NULL; }
	sincos_2pibyn(2 * plan->n, tmp);
	plan->bk[0] = 1;
	plan->bk[1] = 0;

	size_t coeff = 0;
	for (size_t m = 1; m < plan->n; ++m)
	{
		coeff += 2 * m - 1;
		if (coeff >= 2 * plan->n) coeff -= 2 * plan->n;
		plan->bk[2 * m] = tmp[2 * coeff];
		plan->bk[2 * m + 1] = tmp[2 * coeff + 1];
	}

	/* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
	double xn2 = 1. / plan->n2;
	plan->bkf[0] = plan->bk[0] * xn2;
	plan->bkf[1] = plan->bk[1] * xn2;
	for (size_t m = 2; m < 2 * plan->n; m += 2)
	{
		plan->bkf[m] = plan->bkf[2 * plan->n2 - m] = plan->bk[m] * xn2;
		plan->bkf[m + 1] = plan->bkf[2 * plan->n2 - m + 1] = plan->bk[m + 1] * xn2;
	}
	for (size_t m = 2 * plan->n; m <= (2 * plan->n2 - 2 * plan->n + 1); ++m)
		plan->bkf[m] = 0.;
	plan->plan = make_cfftp_plan(plan->n2);
	if (!plan->plan)
	{
		DEALLOC(tmp); DEALLOC(plan->mem); DEALLOC(plan); return NULL;
	}
	if (cfftp_forward(plan->plan, plan->bkf, 1.) != 0)
	{
		DEALLOC(tmp); DEALLOC(plan->mem); DEALLOC(plan); return NULL;
	}
	DEALLOC(tmp);

	return plan;
}

NOINLINE static void destroy_fftblue_plan(fftblue_plan plan)
{
	DEALLOC(plan->mem);
	destroy_cfftp_plan(plan->plan);
	DEALLOC(plan);
}

NOINLINE WARN_UNUSED_RESULT
static int fftblue_fft(fftblue_plan plan, double c[], int isign, double fct)
{
	size_t n = plan->n;
	size_t n2 = plan->n2;
	double* bk = plan->bk;
	double* bkf = plan->bkf;
	double* akf = RALLOC(double, 2 * n2);
	if (!akf) return -1;

	/* initialize a_k and FFT it */
	if (isign > 0)
		for (size_t m = 0; m < 2 * n; m += 2)
		{
			akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
			akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
		}
	else
		for (size_t m = 0; m < 2 * n; m += 2)
		{
			akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
			akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
		}
	for (size_t m = 2 * n; m < 2 * n2; ++m)
		akf[m] = 0;

	if (cfftp_forward(plan->plan, akf, fct) != 0)
	{
		DEALLOC(akf); return -1;
	}

	/* do the convolution */
	if (isign > 0)
		for (size_t m = 0; m < 2 * n2; m += 2)
		{
			double im = -akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
			akf[m] = akf[m] * bkf[m] + akf[m + 1] * bkf[m + 1];
			akf[m + 1] = im;
		}
	else
		for (size_t m = 0; m < 2 * n2; m += 2)
		{
			double im = akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
			akf[m] = akf[m] * bkf[m] - akf[m + 1] * bkf[m + 1];
			akf[m + 1] = im;
		}

	/* inverse FFT */
	if (cfftp_backward(plan->plan, akf, 1.) != 0)
	{
		DEALLOC(akf); return -1;
	}

	/* multiply by b_k */
	if (isign > 0)
		for (size_t m = 0; m < 2 * n; m += 2)
		{
			c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
			c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
		}
	else
		for (size_t m = 0; m < 2 * n; m += 2)
		{
			c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
			c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
		}
	DEALLOC(akf);
	return 0;
}

WARN_UNUSED_RESULT
static int cfftblue_backward(fftblue_plan plan, double c[], double fct)
{
	return fftblue_fft(plan, c, 1, fct);
}

WARN_UNUSED_RESULT
static int cfftblue_forward(fftblue_plan plan, double c[], double fct)
{
	return fftblue_fft(plan, c, -1, fct);
}

WARN_UNUSED_RESULT
static int rfftblue_backward(fftblue_plan plan, double c[], double fct)
{
	size_t n = plan->n;
	double* tmp = RALLOC(double, 2 * n);
	if (!tmp) return -1;
	tmp[0] = c[0];
	tmp[1] = 0.;
	memcpy(tmp + 2, c + 1, (n - 1) * sizeof(double));
	if ((n & 1) == 0) tmp[n + 1] = 0.;
	for (size_t m = 2; m < n; m += 2)
	{
		tmp[2 * n - m] = tmp[m];
		tmp[2 * n - m + 1] = -tmp[m + 1];
	}
	if (fftblue_fft(plan, tmp, 1, fct) != 0)
	{
		DEALLOC(tmp); return -1;
	}
	for (size_t m = 0; m < n; ++m)
		c[m] = tmp[2 * m];
	DEALLOC(tmp);
	return 0;
}

WARN_UNUSED_RESULT
static int rfftblue_forward(fftblue_plan plan, double c[], double fct)
{
	size_t n = plan->n;
	double* tmp = RALLOC(double, 2 * n);
	if (!tmp) return -1;
	for (size_t m = 0; m < n; ++m)
	{
		tmp[2 * m] = c[m];
		tmp[2 * m + 1] = 0.;
	}

	if (fftblue_fft(plan, tmp, -1, fct) != 0)
	{
		DEALLOC(tmp); return -1;
	}
	c[0] = tmp[0];
	memcpy(c + 1, tmp + 2, (n - 1) * sizeof(double));
	DEALLOC(tmp);
	return 0;
}

typedef struct cfft_plan_i
{
	cfftp_plan packplan;
	fftblue_plan blueplan;
} cfft_plan_i;

cfft_plan make_cfft_plan(size_t length)
{
	if (length == 0) return NULL;
	cfft_plan plan = RALLOC(cfft_plan_i, 1);
	if (!plan) return NULL;
	plan->blueplan = 0;
	plan->packplan = 0;
	if ((length < 50) || (largest_prime_factor(length) <= sqrt(length)))
	{
		plan->packplan = make_cfftp_plan(length);
		if (!plan->packplan) { DEALLOC(plan); return NULL; }
		return plan;
	}
	double comp1 = cost_guess(length);
	double comp2 = 2 * cost_guess(good_size(2 * length - 1));
	comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
	if (comp2 < comp1) // use Bluestein
	{
		plan->blueplan = make_fftblue_plan(length);
		if (!plan->blueplan) { DEALLOC(plan); return NULL; }
	}
	else
	{
		plan->packplan = make_cfftp_plan(length);
		if (!plan->packplan) { DEALLOC(plan); return NULL; }
	}
	return plan;
}

void destroy_cfft_plan(cfft_plan plan)
{
	if (plan->blueplan)
		destroy_fftblue_plan(plan->blueplan);
	if (plan->packplan)
		destroy_cfftp_plan(plan->packplan);
	DEALLOC(plan);
}

WARN_UNUSED_RESULT int cfft_backward(cfft_plan plan, double c[], double fct)
{
	if (plan->packplan)
		return cfftp_backward(plan->packplan, c, fct);
	// if (plan->blueplan)
	return cfftblue_backward(plan->blueplan, c, fct);
}

WARN_UNUSED_RESULT int cfft_forward(cfft_plan plan, double c[], double fct)
{
	if (plan->packplan)
		return cfftp_forward(plan->packplan, c, fct);
	// if (plan->blueplan)
	return cfftblue_forward(plan->blueplan, c, fct);
}

typedef struct rfft_plan_i
{
	rfftp_plan packplan;
	fftblue_plan blueplan;
} rfft_plan_i;

rfft_plan make_rfft_plan(size_t length)
{
	if (length == 0) return NULL;
	rfft_plan plan = RALLOC(rfft_plan_i, 1);
	if (!plan) return NULL;
	plan->blueplan = 0;
	plan->packplan = 0;
	if ((length < 50) || (largest_prime_factor(length) <= sqrt(length)))
	{
		plan->packplan = make_rfftp_plan(length);
		if (!plan->packplan) { DEALLOC(plan); return NULL; }
		return plan;
	}
	double comp1 = 0.5 * cost_guess(length);
	double comp2 = 2 * cost_guess(good_size(2 * length - 1));
	comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
	if (comp2 < comp1) // use Bluestein
	{
		plan->blueplan = make_fftblue_plan(length);
		if (!plan->blueplan) { DEALLOC(plan); return NULL; }
	}
	else
	{
		plan->packplan = make_rfftp_plan(length);
		if (!plan->packplan) { DEALLOC(plan); return NULL; }
	}
	return plan;
}

void destroy_rfft_plan(rfft_plan plan)
{
	if (plan->blueplan)
		destroy_fftblue_plan(plan->blueplan);
	if (plan->packplan)
		destroy_rfftp_plan(plan->packplan);
	DEALLOC(plan);
}

size_t rfft_length(rfft_plan plan)
{
	if (plan->packplan) return plan->packplan->length;
	return plan->blueplan->n;
}

size_t cfft_length(cfft_plan plan)
{
	if (plan->packplan) return plan->packplan->length;
	return plan->blueplan->n;
}

WARN_UNUSED_RESULT int rfft_backward(rfft_plan plan, double c[], double fct)
{
	if (plan->packplan)
		return rfftp_backward(plan->packplan, c, fct);
	else // if (plan->blueplan)
		return rfftblue_backward(plan->blueplan, c, fct);
}

WARN_UNUSED_RESULT int rfft_forward(rfft_plan plan, double c[], double fct)
{
	if (plan->packplan)
		return rfftp_forward(plan->packplan, c, fct);
	else // if (plan->blueplan)
		return rfftblue_forward(plan->blueplan, c, fct);
}
