using System.Runtime.InteropServices.ComTypes;

namespace PocketFFT
{
    using System;
    using System.Buffers.Binary;
    using System.Collections;
    using System.Collections.Generic;
    using System.Runtime.InteropServices;

    // https://gitlab.mpcdf.mpg.de/mtr/pocketfft

    // compare also https://www.shoup.net/PGFFT/
    // https://bitbucket.org/jpommier/pffft/src/master/
    public static class PocketFFTBaseline
    {
        private const bool NAIL_TEST = false;

        private static void NailTestPrintComplexArray(ReadOnlySpan<cmplx> val)
        {
            for (int c = 0; c < val.Length; c++)
            {
                if (NAIL_TEST) NailTestPrintComplex(val[c]);
            }
        }

        private static void NailTestPrintComplex(cmplx val)
        {
            if (NAIL_TEST) Console.WriteLine(string.Format("0x{0:X16}r|0x{1:X16}i",
                new DoubleLayout(val.r).LVal,
                new DoubleLayout(val.i).LVal).ToLowerInvariant());
        }

        private static void NailTestPrintDoubleArray(ReadOnlySpan<double> val)
        {
            int col = 0;
            for (int c = 0; c < val.Length; c++)
            {
                if (NAIL_TEST) NailTestPrintDouble(val[c]);
                if (c != (val.Length - 1))
                {
                    Console.Write(", ");
                }
                if (++col > 4)
                {
                    if (NAIL_TEST) Console.WriteLine();
                    col = 0;
                }
            }

            if (NAIL_TEST) Console.WriteLine();
        }

        private static void NailTestPrintDouble(double val)
        {
            Console.Write(string.Format("0x{0:X}", new DoubleLayout(val).LVal));
        }

        [StructLayout(LayoutKind.Explicit)]
        private struct DoubleLayout
        {
            [FieldOffset(0)]
            public double DVal;

            [FieldOffset(0)]
            public long LVal;

            public DoubleLayout(double d)
            {
                LVal = 0;
                DVal = d;
            }
        }


        //#define SWAP(a,b,type)
        //    do { type tmp_=(a); (a)=(b); (b)=tmp_; } while(0)
        private static void Swap<T>(ref T a, ref T b)
        {
            T tmp = a;
            a = b;
            b = tmp;
        }

        private static double fma(double a, double b, double c)
        {
#if NET8_0_OR_GREATER
            return Math.FusedMultiplyAdd(a, b, c);
#else
            return (a * b) + c;
#endif
        }

        // adapted from https://stackoverflow.com/questions/42792939/
        // CAUTION: this function only works for arguments in the range [-0.25; 0.25]!
        private static void my_sincosm1pi(double a, Span<double> res)
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

        private static void calc_first_octant(int den, Span<double> res)
        {
            Span<double> cs = stackalloc double[2];
            int n = (den + 4) >> 3;
            if (n == 0)
            {
                return;
            }

            res[0] = 1.0;
            res[1] = 0.0;
            if (n == 1)
            {
                return;
            }

            int l1 = (int)Math.Sqrt(n);
            for (int i = 1; i < l1; ++i)
            {
                my_sincosm1pi((2.0 * i) / den, res.Slice(2 * i));
            }

            int start = l1;
            while (start < n)
            {
                my_sincosm1pi((2.0 * start) / den, cs);
                res[2 * start] = cs[0] + 1.0;
                res[2 * start + 1] = cs[1];

                int end = l1;
                if (start + end > n)
                {
                    end = n - start;
                }

                for (int i = 1; i < end; ++i)
                {
                    double csx0 = res[2 * i];
                    double csx1 = res[2 * i + 1];
                    res[2 * (start + i)] = ((cs[0] * csx0 - cs[1] * csx1 + cs[0]) + csx0) + 1.0;
                    res[2 * (start + i) + 1] = (cs[0] * csx1 + cs[1] * csx0) + cs[1] + csx1;
                }

                start += l1;
            }

            for (int i = 1; i < l1; ++i)
            {
                res[2 * i] += 1.0;
            }
        }

        private static void calc_first_quadrant(int n, Span<double> res)
        {
            Span<double> p = res.Slice(n);
            calc_first_octant(n << 1, p);
            int ndone = (n + 2) >> 2;
            int i = 0, idx1 = 0, idx2 = 2 * ndone - 2;
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

        private static void calc_first_half(int n, Span<double> res)
        {
            int ndone = (n + 1) >> 1;
            Span<double> p = res.Slice(n - 1);
            calc_first_octant(n << 2, p);
            int i4 = 0, inp = n, i = 0;

            for (; i4 <= inp - i4; ++i, i4 += 4) // octant 0
            {
                res[2 * i] = p[2 * i4]; res[2 * i + 1] = p[2 * i4 + 1];
            }

            for (; i4 - inp <= 0; ++i, i4 += 4) // octant 1
            {
                int xm = inp - i4;
                res[2 * i] = p[2 * xm + 1]; res[2 * i + 1] = p[2 * xm];
            }

            for (; i4 <= 3 * inp - i4; ++i, i4 += 4) // octant 2
            {
                int xm = i4 - inp;
                res[2 * i] = -p[2 * xm + 1]; res[2 * i + 1] = p[2 * xm];
            }

            for (; i < ndone; ++i, i4 += 4) // octant 3
            {
                int xm = 2 * inp - i4;
                res[2 * i] = -p[2 * xm]; res[2 * i + 1] = p[2 * xm + 1];
            }
        }

        private static void fill_first_quadrant(int n, Span<double> res)
        {
            const double hsqt2 = 0.707106781186547524400844362104849;
            int quart = n >> 2;

            if ((n & 7) == 0)
            {
                res[quart] = res[quart + 1] = hsqt2;
            }

            for (int i = 2, j = 2 * quart - 2; i < quart; i += 2, j -= 2)
            {
                res[j] = res[i + 1];
                res[j + 1] = res[i];
            }
        }

        private static void fill_first_half(int n, Span<double> res)
        {
            int half = n >> 1;
            if ((n & 3) == 0)
            {
                for (int i = 0; i < half; i += 2)
                {
                    res[i + half] = -res[i + 1];
                    res[i + half + 1] = res[i];
                }
            }
            else
            {
                for (int i = 2, j = 2 * half - 2; i < half; i += 2, j -= 2)
                {
                    res[j] = -res[i];
                    res[j + 1] = res[i + 1];
                }
            }
        }

        private static void fill_second_half(int n, Span<double> res)
        {
            if ((n & 1) == 0)
            {
                for (int i = 0; i < n; ++i)
                {
                    res[i + n] = -res[i];
                }
            }
            else
            {
                for (int i = 2, j = 2 * n - 2; i < n; i += 2, j -= 2)
                {
                    res[j] = res[i];
                    res[j + 1] = -res[i + 1];
                }
            }
        }

        private static void sincos_2pibyn_half(int n, Span<double> res)
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
            {
                calc_first_half(n, res);
            }
        }

        private static void sincos_2pibyn(int n, Span<double> res)
        {
            sincos_2pibyn_half(n, res);
            fill_second_half(n, res);
        }

        private static int largest_prime_factor(int n)
        {
            int res = 1;
            int tmp;
            while (((tmp = (n >> 1)) << 1) == n)
            {
                res = 2;
                n = tmp;
            }

            int limit = (int)Math.Sqrt(n + 0.01);
            for (int x = 3; x <= limit; x += 2)
            {
                while (((tmp = (n / x)) * x) == n)
                {
                    res = x;
                    n = tmp;
                    limit = (int)Math.Sqrt(n + 0.01);
                }
            }

            if (n > 1)
            {
                res = n;
            }

            return res;
        }

        private static double cost_guess(int n)
        {
            const double lfp = 1.1; // penalty for non-hardcoded larger factors
            int ni = n;
            double result = 0.0;
            int tmp;
            while (((tmp = (n >> 1)) << 1) == n)
            {
                result += 2;
                n = tmp;
            }

            int limit = (int)Math.Sqrt(n + 0.01);
            for (int x = 3; x <= limit; x += 2)
            {
                while ((tmp = (n / x)) * x == n)
                {
                    result += (x <= 5) ? x : lfp * x; // penalize larger prime factors
                    n = tmp;
                    limit = (int)Math.Sqrt(n + 0.01);
                }
            }

            if (n > 1)
            {
                result += (n <= 5) ? n : lfp * n;
            }

            return result * ni;
        }

        /* returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n */
        private static int good_size(int n)
        {
            if (n <= 6) return n;

            int bestfac = 2 * n;
            for (int f2 = 1; f2 < bestfac; f2 *= 2)
                for (int f23 = f2; f23 < bestfac; f23 *= 3)
                    for (int f235 = f23; f235 < bestfac; f235 *= 5)
                        for (int f2357 = f235; f2357 < bestfac; f2357 *= 7)
                            for (int f235711 = f2357; f235711 < bestfac; f235711 *= 11)
                                if (f235711 >= n) bestfac = f235711;
            return bestfac;
        }

        public struct cmplx
        {
            public double r;
            public double i;

            public cmplx(double r, double i)
            {
                this.r = r;
                this.i = i;
            }
        }

        public struct cfftp_fctdata
        {
            public int fct;
            public int tw; // originally these were cmplx* pointers, but they just indexed into plan.mem
            public int tws; // so they have been replaced by basic integer indexes

            public cfftp_fctdata(int fct, int tw, int tws)
            {
                this.fct = fct;
                this.tw = tw;
                this.tws = tws;
            }
        }

        private const int NFCT = 25;

#if NET9_0_OR_GREATER
        [System.Runtime.CompilerServices.InlineArray(NFCT)]
        public struct cfftp_fctdata_array
        {
            public cfftp_fctdata data;
        }

        public class cfftp_plan
        {
            public int length;
            public int nfct;
            public cmplx[] mem;
            public cfftp_fctdata_array fct;
        }
#else
        public class cfftp_plan
        {
            public int length;
            public int nfct;
            public cmplx[] mem;
            public cfftp_fctdata[] fct; // [NFCT]

            public IEnumerable<int> Factors()
            {
                foreach (var data in fct)
                {
                    if (data.fct == 0)
                    {
                        yield break;
                    }

                    yield return data.fct;
                }
            }
        }
#endif

        // #define PMC(a,b,c,d) { a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
        private static void PMC(ref cmplx a, ref cmplx b, ref cmplx c, ref cmplx d)
        {
            a.r = c.r + d.r;
            a.i = c.i + d.i;
            b.r = c.r - d.r;
            b.i = c.i - d.i;
        }

        //#define ADDC(a,b,c) { a.r=b.r+c.r; a.i=b.i+c.i; }
        private static void ADDC(ref cmplx a, ref cmplx b, ref cmplx c)
        {
            a.r = b.r + c.r;
            a.i = b.i + c.i;
        }

        //#define SCALEC(a,b) { a.r*=b; a.i*=b; }
        private static void SCALEC(ref cmplx a, double b)
        {
            a.r *= b;
            a.i *= b;
        }

        //#define ROT90(a) { double tmp_=a.r; a.r=-a.i; a.i=tmp_; }
        private static void ROT90(ref cmplx a)
        {
            double tmp = a.r;
            a.r = -a.i;
            a.i = tmp;
        }

        //#define ROTM90(a) { double tmp_=-a.r; a.r=a.i; a.i=tmp_; }
        private static void ROTM90(ref cmplx a)
        {
            double tmp = -a.r;
            a.r = a.i;
            a.i = tmp;
        }

        //#define A_EQ_B_MUL_C(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
        private static void A_EQ_B_MUL_C(ref cmplx a, ref cmplx b, ref cmplx c)
        {
            /* a = b*c */
            a.r = b.r * c.r - b.i * c.i;
            a.i = b.r * c.i + b.i * c.r;
        }

        //#define A_EQ_CB_MUL_C(a,b,c) { a.r=b.r*c.r+b.i*c.i; a.i=b.r*c.i-b.i*c.r; }
        private static void A_EQ_CB_MUL_C(ref cmplx a, ref cmplx b, ref cmplx c)
        {
            /* a = conj(b)*c*/
            a.r = b.r * c.r + b.i * c.i;
            a.i = b.r * c.i - b.i * c.r;
        }

        //#define PMSIGNC(a,b,c,d) { a.r=c.r+sign*d.r; a.i=c.i+sign*d.i; b.r=c.r-sign*d.r; b.i=c.i-sign*d.i; }
        private static void PMSIGNC(ref cmplx a, ref cmplx b, ref cmplx c, ref cmplx d, int sign)
        {
            a.r = c.r + sign * d.r;
            a.i = c.i + sign * d.i;
            b.r = c.r - sign * d.r;
            b.i = c.i - sign * d.i;
        }

        //#define MULPMSIGNC(a,b,c) { a.r=b.r*c.r-sign*b.i*c.i; a.i=b.r*c.i+sign*b.i*c.r; }
        private static void MULPMSIGNC(ref cmplx a, ref cmplx b, ref cmplx c, int sign)
        {
            /* a = b*c */
            a.r = b.r * c.r - sign * b.i * c.i;
            a.i = b.r * c.i + sign * b.i * c.r;
        }

        //#define MULPMSIGNCEQ(a,b) { double xtmp=a.r; a.r=b.r*a.r-sign*b.i*a.i; a.i=b.r*a.i+sign*b.i*xtmp; }
        private static void MULPMSIGNCEQ(ref cmplx a, ref cmplx b, int sign)
        {
            /* a *= b */
            double xtmp = a.r;
            a.r = b.r * a.r - sign * b.i * a.i;
            a.i = b.r * a.i + sign * b.i * xtmp;
        }

        private static void pass2b_opt(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;

            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    PMC(ref ch[(0) + ido * ((k) + l1 * (0))],
                        ref ch[(0) + ido * ((k) + l1 * (1))],
                        ref cc[(0) + ido * ((0) + cdim * (k))],
                        ref cc[(0) + ido * ((1) + cdim * (k))]);
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    PMC(ref ch[(0) + ido * ((k) + l1 * (0))],
                        ref ch[(0) + ido * ((k) + l1 * (1))],
                        ref cc[(0) + ido * ((0) + cdim * (k))],
                        ref cc[(0) + ido * ((1) + cdim * (k))]);

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t = default;
                        PMC(ref ch[(i) + ido * ((k) + l1 * (0))],
                            ref t,
                            ref cc[(i) + ido * ((0) + cdim * (k))],
                            ref cc[(i) + ido * ((1) + cdim * (k))]);
                        A_EQ_B_MUL_C(
                            ref ch[(i) + ido * ((k) + l1 * (1))],
                            ref wa[(i) - 1 + (0) * (ido - 1)],
                            ref t);
                    }
                }
            }
        }

        private static void pass2b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;
            if (NAIL_TEST) Console.WriteLine("start pass2b");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t;
                        ch[(i) + ido * ((k) + l1 * (0))].r = cc[(i) + ido * ((0) + cdim * (k))].r + cc[(i) + ido * ((1) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = cc[(i) + ido * ((0) + cdim * (k))].i + cc[(i) + ido * ((1) + cdim * (k))].i;
                        t.r = cc[(i) + ido * ((0) + cdim * (k))].r - cc[(i) + ido * ((1) + cdim * (k))].r;
                        t.i = cc[(i) + ido * ((0) + cdim * (k))].i - cc[(i) + ido * ((1) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (0) * (ido - 1)].r * t.r - wa[(i) - 1 + (0) * (ido - 1)].i * t.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (0) * (ido - 1)].r * t.i + wa[(i) - 1 + (0) * (ido - 1)].i * t.r;
                    }
                }
            }
        }

        private static void pass2f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 2;
            if (NAIL_TEST) Console.WriteLine("start pass2f");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))].r = cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i;
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t;
                        ch[(i) + ido * ((k) + l1 * (0))].r = cc[(i) + ido * ((0) + cdim * (k))].r + cc[(i) + ido * ((1) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = cc[(i) + ido * ((0) + cdim * (k))].i + cc[(i) + ido * ((1) + cdim * (k))].i;
                        t.r = cc[(i) + ido * ((0) + cdim * (k))].r - cc[(i) + ido * ((1) + cdim * (k))].r;
                        t.i = cc[(i) + ido * ((0) + cdim * (k))].i - cc[(i) + ido * ((1) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (0) * (ido - 1)].r * t.r + wa[(i) - 1 + (0) * (ido - 1)].i * t.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (0) * (ido - 1)].r * t.i - wa[(i) - 1 + (0) * (ido - 1)].i * t.r;
                    }
                }
            }
        }

        private static void pass3b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 3;
            const double tw1r = -0.5, tw1i = 0.86602540378443864676;
            if (NAIL_TEST) Console.WriteLine("start pass3b");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r; ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                        t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                        t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                        ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                        cmplx ca, cb;
                        ca.r = t0.r + tw1r * t1.r;
                        ca.i = t0.i + tw1r * t1.i;
                        cb.i = tw1i * t2.r;
                        cb.r = -(tw1i * t2.i);
                        ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                        ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                        ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                        ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                    }

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2;
                        t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((2) + cdim * (k))].r;
                        t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((2) + cdim * (k))].i;
                        t2.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((2) + cdim * (k))].r;
                        t2.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((2) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                        cmplx ca, cb, da, db;
                        ca.r = t0.r + tw1r * t1.r;
                        ca.i = t0.i + tw1r * t1.i;
                        cb.i = tw1i * t2.r;
                        cb.r = -(tw1i * t2.i);
                        da.r = ca.r + cb.r;
                        da.i = ca.i + cb.i;
                        db.r = ca.r - cb.r;
                        db.i = ca.i - cb.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }

        private static void pass3f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 3;
            const double tw1r = -0.5, tw1i = -0.86602540378443864676;

            if (NAIL_TEST) Console.WriteLine("start pass3f");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2;
                    t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((2) + cdim * (k))].r;
                    t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((2) + cdim * (k))].i;
                    t2.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((2) + cdim * (k))].r;
                    t2.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((2) + cdim * (k))].i;
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i;
                    cmplx ca, cb;
                    ca.r = t0.r + tw1r * t1.r;
                    ca.i = t0.i + tw1r * t1.i;
                    cb.i = tw1i * t2.r;
                    cb.r = -(tw1i * t2.i);
                    ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = ca.r - cb.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = ca.i - cb.i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t3 = cc[(i) + ido * ((0) + cdim * (k))], t4, t5;
                        t4.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((2) + cdim * (k))].r;
                        t4.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((2) + cdim * (k))].i;
                        t5.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((2) + cdim * (k))].r;
                        t5.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((2) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (0))].r = t3.r + t4.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t3.i + t4.i;
                        cmplx ca2, cb2, da, db;
                        ca2.r = t3.r + tw1r * t4.r;
                        ca2.i = t3.i + tw1r * t4.i;
                        cb2.i = tw1i * t5.r;
                        cb2.r = -(tw1i * t5.i);
                        da.r = ca2.r + cb2.r;
                        da.i = ca2.i + cb2.i;
                        db.r = ca2.r - cb2.r;
                        db.i = ca2.i - cb2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * db.r;
                    }
                }
            }
        }

        private static void pass4b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 4;
            if (NAIL_TEST) Console.WriteLine("start pass4b");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
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
                    ROT90(ref t4);
        
            ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
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
                    ROT90(ref t4);
        
                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;

                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx c2, c3, c4, t5, t6, t7, t8;
                        cmplx cc0 = cc[(i) + ido * ((0) + cdim * (k))],
                            cc1 = cc[(i) + ido * ((1) + cdim * (k))],
                            cc2 = cc[(i) + ido * ((2) + cdim * (k))],
                            cc3 = cc[(i) + ido * ((3) + cdim * (k))];
                        t6.r = cc0.r + cc2.r;
                        t6.i = cc0.i + cc2.i;
                        t5.r = cc0.r - cc2.r;
                        t5.i = cc0.i - cc2.i;
                        t7.r = cc1.r + cc3.r;
                        t7.i = cc1.i + cc3.i;
                        t8.r = cc1.r - cc3.r;
                        t8.i = cc1.i - cc3.i;
                        ROT90(ref t8);

                        cmplx wa0 = wa[(i) - 1 + (0) * (ido - 1)],
                        wa1 = wa[(i) - 1 + (1) * (ido - 1)],
                        wa2 = wa[(i) - 1 + (2) * (ido - 1)];
                        ch[(i) + ido * ((k) + l1 * (0))].r = t6.r + t7.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t6.i + t7.i;
                        c3.r = t6.r - t7.r;
                        c3.i = t6.i - t7.i;
                        c2.r = t5.r + t8.r;
                        c2.i = t5.i + t8.i;
                        c4.r = t5.r - t8.r;
                        c4.i = t5.i - t8.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa0.r * c2.r - wa0.i * c2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa0.r * c2.i + wa0.i * c2.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa1.r * c3.r - wa1.i * c3.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa1.r * c3.i + wa1.i * c3.r;
                        ch[(i) + ido * ((k) + l1 * (3))].r = wa2.r * c4.r - wa2.i * c4.i;
                        ch[(i) + ido * ((k) + l1 * (3))].i = wa2.r * c4.i + wa2.i * c4.r;
                    }
                }
            }
        }

        private static void pass4f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 4;

            if (NAIL_TEST) Console.WriteLine("start pass4f");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
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
                    ROTM90(ref t4);
        
                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
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
                    ROTM90(ref t4);
        
                    ch[(0) + ido * ((k) + l1 * (0))].r = t2.r + t3.r;
                    ch[(0) + ido * ((k) + l1 * (0))].i = t2.i + t3.i;
                    ch[(0) + ido * ((k) + l1 * (2))].r = t2.r - t3.r;
                    ch[(0) + ido * ((k) + l1 * (2))].i = t2.i - t3.i;
                    ch[(0) + ido * ((k) + l1 * (1))].r = t1.r + t4.r;
                    ch[(0) + ido * ((k) + l1 * (1))].i = t1.i + t4.i;
                    ch[(0) + ido * ((k) + l1 * (3))].r = t1.r - t4.r;
                    ch[(0) + ido * ((k) + l1 * (3))].i = t1.i - t4.i;
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx c2, c3, c4, t5, t6, t7, t8;
                        cmplx cc0 = cc[(i) + ido * ((0) + cdim * (k))],
                            cc1 = cc[(i) + ido * ((1) + cdim * (k))],
                            cc2 = cc[(i) + ido * ((2) + cdim * (k))],
                            cc3 = cc[(i) + ido * ((3) + cdim * (k))];
                        t6.r = cc0.r + cc2.r;
                        t6.i = cc0.i + cc2.i;
                        t5.r = cc0.r - cc2.r;
                        t5.i = cc0.i - cc2.i;
                        t7.r = cc1.r + cc3.r;
                        t7.i = cc1.i + cc3.i;
                        t8.r = cc1.r - cc3.r;
                        t8.i = cc1.i - cc3.i;
                        ROTM90(ref t8);

                        cmplx wa0 = wa[(i) - 1 + (0) * (ido - 1)],
                        wa1 = wa[(i) - 1 + (1) * (ido - 1)],
                        wa2 = wa[(i) - 1 + (2) * (ido - 1)];
                        ch[(i) + ido * ((k) + l1 * (0))].r = t6.r + t7.r;
                        ch[(i) + ido * ((k) + l1 * (0))].i = t6.i + t7.i;
                        c3.r = t6.r - t7.r;
                        c3.i = t6.i - t7.i;
                        c2.r = t5.r + t8.r;
                        c2.i = t5.i + t8.i;
                        c4.r = t5.r - t8.r;
                        c4.i = t5.i - t8.i;
                        ch[(i) + ido * ((k) + l1 * (1))].r = wa0.r * c2.r + wa0.i * c2.i;
                        ch[(i) + ido * ((k) + l1 * (1))].i = wa0.r * c2.i - wa0.i * c2.r;
                        ch[(i) + ido * ((k) + l1 * (2))].r = wa1.r * c3.r + wa1.i * c3.i;
                        ch[(i) + ido * ((k) + l1 * (2))].i = wa1.r * c3.i - wa1.i * c3.r;
                        ch[(i) + ido * ((k) + l1 * (3))].r = wa2.r * c4.r + wa2.i * c4.i;
                        ch[(i) + ido * ((k) + l1 * (3))].i = wa2.r * c4.i - wa2.i * c4.r;
                    }
                }
            }
        }

        static void pass5b(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {

            const int cdim = 5;
            const double tw1r = 0.3090169943749474241,
                tw1i = 0.95105651629515357212,
                tw2r = -0.8090169943749474241,
                tw2i = 0.58778525229247312917;
            if (NAIL_TEST) Console.WriteLine("start pass5b");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));

            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                    {
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    {
                        t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                    {
                        cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                        }
                    }
                }

            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t4.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((3) + cdim * (k))].r; t2.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((3) + cdim * (k))].i; t3.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((3) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.r - wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.i + wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass5f(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa)
        {
            const int cdim = 5;
            const double tw1r = 0.3090169943749474241,
                tw1i = -0.95105651629515357212,
                tw2r = -0.8090169943749474241,
                tw2i = -0.58778525229247312917;

            if (NAIL_TEST) Console.WriteLine("start pass5f");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                    {
                        t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    {
                        t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                    {
                        cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t0 = cc[(0) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t1.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t4.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((3) + cdim * (k))].r; t2.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((3) + cdim * (k))].i; t3.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((3) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(0) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (3))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t0 = cc[(i) + ido * ((0) + cdim * (k))], t1, t2, t3, t4;
                        {
                            t1.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t1.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t4.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        {
                            t2.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((3) + cdim * (k))].r; t2.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((3) + cdim * (k))].i; t3.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((3) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((3) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t0.r + t1.r + t2.r; ch[(i) + ido * ((k) + l1 * (0))].i = t0.i + t1.i + t2.i;
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw1r * t1.r + tw2r * t2.r; ca.i = t0.i + tw1r * t1.i + tw2r * t2.i; cb.i = +tw1i * t4.r + tw2i * t3.r; cb.r = -(+tw1i * t4.i + tw2i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx ca, cb, da, db; ca.r = t0.r + tw2r * t1.r + tw1r * t2.r; ca.i = t0.i + tw2r * t1.i + tw1r * t2.i; cb.i = +tw2i * t4.r - tw1i * t3.r; cb.r = -(+tw2i * t4.i - tw1i * t3.i);
                            {
                                da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r + wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i - wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.r + wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * db.i - wa[(i) - 1 + (3 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass7(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, int sign)
        {
            const int cdim = 7;
            double tw1r = 0.623489801858733530525,
                tw1i = sign * 0.7818314824680298087084,
                tw2r = -0.222520933956314404289,
                tw2i = sign * 0.9749279121818236070181,
                tw3r = -0.9009688679024191262361,
                tw3i = sign * 0.4338837391175581204758;

            if (NAIL_TEST) Console.WriteLine("start pass7");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                    {
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                    }
                    {
                        t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
                    }
                    {
                        t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                    {
                        cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                        {
                            t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((5) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((5) + cdim * (k))].i; t6.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((5) + cdim * (k))].r; t6.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((5) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((4) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((4) + cdim * (k))].i; t5.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((4) + cdim * (k))].r; t5.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((4) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                        {
                            cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (5))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (4))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t1 = cc[(i) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7;
                        {
                            t2.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((6) + cdim * (k))].r; t2.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((6) + cdim * (k))].i; t7.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((6) + cdim * (k))].r; t7.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((6) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((5) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((5) + cdim * (k))].i; t6.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((5) + cdim * (k))].r; t6.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((5) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(i) + ido * ((3) + cdim * (k))].r + cc[(i) + ido * ((4) + cdim * (k))].r; t4.i = cc[(i) + ido * ((3) + cdim * (k))].i + cc[(i) + ido * ((4) + cdim * (k))].i; t5.r = cc[(i) + ido * ((3) + cdim * (k))].r - cc[(i) + ido * ((4) + cdim * (k))].r; t5.i = cc[(i) + ido * ((3) + cdim * (k))].i - cc[(i) + ido * ((4) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r; ch[(i) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i;
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i; cb.i = +tw1i * t7.r + tw2i * t6.r + tw3i * t5.r; cb.r = -(+tw1i * t7.i + tw2i * t6.i + tw3i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (6))].r = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (6))].i = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw3r * t3.r + tw1r * t4.r; ca.i = t1.i + tw2r * t2.i + tw3r * t3.i + tw1r * t4.i; cb.i = +tw2i * t7.r - tw3i * t6.r - tw1i * t5.r; cb.r = -(+tw2i * t7.i - tw3i * t6.i - tw1i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (5))].r = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (5))].i = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw1r * t3.r + tw2r * t4.r; ca.i = t1.i + tw3r * t2.i + tw1r * t3.i + tw2r * t4.i; cb.i = +tw3i * t7.r - tw1i * t6.r + tw2i * t5.r; cb.r = -(+tw3i * t7.i - tw1i * t6.i + tw2i * t5.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void pass11(int ido, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, int sign)
        {
            const int cdim = 11;
            double tw1r = 0.8412535328311811688618,
                tw1i = sign * 0.5406408174555975821076,
                tw2r = 0.4154150130018864255293,
                tw2i = sign * 0.9096319953545183714117,
                tw3r = -0.1423148382732851404438,
                tw3i = sign * 0.9898214418809327323761,
                tw4r = -0.6548607339452850640569,
                tw4i = sign * 0.755749574354258283774,
                tw5r = -0.9594929736144973898904,
                tw5i = sign * 0.2817325568414296977114;

            if (NAIL_TEST) Console.WriteLine("start pass11");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            if (ido == 1)
                for (int k = 0; k < l1; ++k)
                {
                    cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                    {
                        t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
                    }
                    {
                        t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
                    }
                    {
                        t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
                    }
                    {
                        t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
                    }
                    {
                        t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                    }
                    ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                    {
                        cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (10))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (9))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (8))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (7))].i = ca.i - cb.i;
                        }
                    }
                    {
                        cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                        {
                            ch[(0) + ido * ((k) + l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                        }
                    }
                }
            else
                for (int k = 0; k < l1; ++k)
                {
                    {
                        cmplx t1 = cc[(0) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                        {
                            t2.r = cc[(0) + ido * ((1) + cdim * (k))].r + cc[(0) + ido * ((10) + cdim * (k))].r; t2.i = cc[(0) + ido * ((1) + cdim * (k))].i + cc[(0) + ido * ((10) + cdim * (k))].i; t11.r = cc[(0) + ido * ((1) + cdim * (k))].r - cc[(0) + ido * ((10) + cdim * (k))].r; t11.i = cc[(0) + ido * ((1) + cdim * (k))].i - cc[(0) + ido * ((10) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(0) + ido * ((2) + cdim * (k))].r + cc[(0) + ido * ((9) + cdim * (k))].r; t3.i = cc[(0) + ido * ((2) + cdim * (k))].i + cc[(0) + ido * ((9) + cdim * (k))].i; t10.r = cc[(0) + ido * ((2) + cdim * (k))].r - cc[(0) + ido * ((9) + cdim * (k))].r; t10.i = cc[(0) + ido * ((2) + cdim * (k))].i - cc[(0) + ido * ((9) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(0) + ido * ((3) + cdim * (k))].r + cc[(0) + ido * ((8) + cdim * (k))].r; t4.i = cc[(0) + ido * ((3) + cdim * (k))].i + cc[(0) + ido * ((8) + cdim * (k))].i; t9.r = cc[(0) + ido * ((3) + cdim * (k))].r - cc[(0) + ido * ((8) + cdim * (k))].r; t9.i = cc[(0) + ido * ((3) + cdim * (k))].i - cc[(0) + ido * ((8) + cdim * (k))].i;
                        }
                        {
                            t5.r = cc[(0) + ido * ((4) + cdim * (k))].r + cc[(0) + ido * ((7) + cdim * (k))].r; t5.i = cc[(0) + ido * ((4) + cdim * (k))].i + cc[(0) + ido * ((7) + cdim * (k))].i; t8.r = cc[(0) + ido * ((4) + cdim * (k))].r - cc[(0) + ido * ((7) + cdim * (k))].r; t8.i = cc[(0) + ido * ((4) + cdim * (k))].i - cc[(0) + ido * ((7) + cdim * (k))].i;
                        }
                        {
                            t6.r = cc[(0) + ido * ((5) + cdim * (k))].r + cc[(0) + ido * ((6) + cdim * (k))].r; t6.i = cc[(0) + ido * ((5) + cdim * (k))].i + cc[(0) + ido * ((6) + cdim * (k))].i; t7.r = cc[(0) + ido * ((5) + cdim * (k))].r - cc[(0) + ido * ((6) + cdim * (k))].r; t7.i = cc[(0) + ido * ((5) + cdim * (k))].i - cc[(0) + ido * ((6) + cdim * (k))].i;
                        }
                        ch[(0) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(0) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                        {
                            cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (1))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (1))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (10))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (10))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (2))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (2))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (9))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (9))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (3))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (3))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (8))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (8))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (4))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (4))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (7))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (7))].i = ca.i - cb.i;
                            }
                        }
                        {
                            cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                            {
                                ch[(0) + ido * ((k) + l1 * (5))].r = ca.r + cb.r; ch[(0) + ido * ((k) + l1 * (5))].i = ca.i + cb.i; ch[(0) + ido * ((k) + l1 * (6))].r = ca.r - cb.r; ch[(0) + ido * ((k) + l1 * (6))].i = ca.i - cb.i;
                            }
                        }
                    }
                    for (int i = 1; i < ido; ++i)
                    {
                        cmplx t1 = cc[(i) + ido * ((0) + cdim * (k))], t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
                        {
                            t2.r = cc[(i) + ido * ((1) + cdim * (k))].r + cc[(i) + ido * ((10) + cdim * (k))].r; t2.i = cc[(i) + ido * ((1) + cdim * (k))].i + cc[(i) + ido * ((10) + cdim * (k))].i; t11.r = cc[(i) + ido * ((1) + cdim * (k))].r - cc[(i) + ido * ((10) + cdim * (k))].r; t11.i = cc[(i) + ido * ((1) + cdim * (k))].i - cc[(i) + ido * ((10) + cdim * (k))].i;
                        }
                        {
                            t3.r = cc[(i) + ido * ((2) + cdim * (k))].r + cc[(i) + ido * ((9) + cdim * (k))].r; t3.i = cc[(i) + ido * ((2) + cdim * (k))].i + cc[(i) + ido * ((9) + cdim * (k))].i; t10.r = cc[(i) + ido * ((2) + cdim * (k))].r - cc[(i) + ido * ((9) + cdim * (k))].r; t10.i = cc[(i) + ido * ((2) + cdim * (k))].i - cc[(i) + ido * ((9) + cdim * (k))].i;
                        }
                        {
                            t4.r = cc[(i) + ido * ((3) + cdim * (k))].r + cc[(i) + ido * ((8) + cdim * (k))].r; t4.i = cc[(i) + ido * ((3) + cdim * (k))].i + cc[(i) + ido * ((8) + cdim * (k))].i; t9.r = cc[(i) + ido * ((3) + cdim * (k))].r - cc[(i) + ido * ((8) + cdim * (k))].r; t9.i = cc[(i) + ido * ((3) + cdim * (k))].i - cc[(i) + ido * ((8) + cdim * (k))].i;
                        }
                        {
                            t5.r = cc[(i) + ido * ((4) + cdim * (k))].r + cc[(i) + ido * ((7) + cdim * (k))].r; t5.i = cc[(i) + ido * ((4) + cdim * (k))].i + cc[(i) + ido * ((7) + cdim * (k))].i; t8.r = cc[(i) + ido * ((4) + cdim * (k))].r - cc[(i) + ido * ((7) + cdim * (k))].r; t8.i = cc[(i) + ido * ((4) + cdim * (k))].i - cc[(i) + ido * ((7) + cdim * (k))].i;
                        }
                        {
                            t6.r = cc[(i) + ido * ((5) + cdim * (k))].r + cc[(i) + ido * ((6) + cdim * (k))].r; t6.i = cc[(i) + ido * ((5) + cdim * (k))].i + cc[(i) + ido * ((6) + cdim * (k))].i; t7.r = cc[(i) + ido * ((5) + cdim * (k))].r - cc[(i) + ido * ((6) + cdim * (k))].r; t7.i = cc[(i) + ido * ((5) + cdim * (k))].i - cc[(i) + ido * ((6) + cdim * (k))].i;
                        }
                        ch[(i) + ido * ((k) + l1 * (0))].r = t1.r + t2.r + t3.r + t4.r + t5.r + t6.r; ch[(i) + ido * ((k) + l1 * (0))].i = t1.i + t2.i + t3.i + t4.i + t5.i + t6.i;
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw1r * t2.r + tw2r * t3.r + tw3r * t4.r + tw4r * t5.r + tw5r * t6.r; ca.i = t1.i + tw1r * t2.i + tw2r * t3.i + tw3r * t4.i + tw4r * t5.i + tw5r * t6.i; cb.i = +tw1i * t11.r + tw2i * t10.r + tw3i * t9.r + tw4i * t8.r + tw5i * t7.r; cb.r = -(+tw1i * t11.i + tw2i * t10.i + tw3i * t9.i + tw4i * t8.i + tw5i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (1))].r = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (1))].i = wa[(i) - 1 + (1 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (1 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (10))].r = wa[(i) - 1 + (10 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (10 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (10))].i = wa[(i) - 1 + (10 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (10 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw2r * t2.r + tw4r * t3.r + tw5r * t4.r + tw3r * t5.r + tw1r * t6.r; ca.i = t1.i + tw2r * t2.i + tw4r * t3.i + tw5r * t4.i + tw3r * t5.i + tw1r * t6.i; cb.i = +tw2i * t11.r + tw4i * t10.r - tw5i * t9.r - tw3i * t8.r - tw1i * t7.r; cb.r = -(+tw2i * t11.i + tw4i * t10.i - tw5i * t9.i - tw3i * t8.i - tw1i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (2))].r = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (2))].i = wa[(i) - 1 + (2 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (2 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (9))].r = wa[(i) - 1 + (9 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (9 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (9))].i = wa[(i) - 1 + (9 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (9 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw3r * t2.r + tw5r * t3.r + tw2r * t4.r + tw1r * t5.r + tw4r * t6.r; ca.i = t1.i + tw3r * t2.i + tw5r * t3.i + tw2r * t4.i + tw1r * t5.i + tw4r * t6.i; cb.i = +tw3i * t11.r - tw5i * t10.r - tw2i * t9.r + tw1i * t8.r + tw4i * t7.r; cb.r = -(+tw3i * t11.i - tw5i * t10.i - tw2i * t9.i + tw1i * t8.i + tw4i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (3))].r = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (3))].i = wa[(i) - 1 + (3 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (3 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (8))].r = wa[(i) - 1 + (8 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (8 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (8))].i = wa[(i) - 1 + (8 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (8 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw4r * t2.r + tw3r * t3.r + tw1r * t4.r + tw5r * t5.r + tw2r * t6.r; ca.i = t1.i + tw4r * t2.i + tw3r * t3.i + tw1r * t4.i + tw5r * t5.i + tw2r * t6.i; cb.i = +tw4i * t11.r - tw3i * t10.r + tw1i * t9.r + tw5i * t8.r - tw2i * t7.r; cb.r = -(+tw4i * t11.i - tw3i * t10.i + tw1i * t9.i + tw5i * t8.i - tw2i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (4))].r = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (4))].i = wa[(i) - 1 + (4 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (4 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (7))].r = wa[(i) - 1 + (7 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (7 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (7))].i = wa[(i) - 1 + (7 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (7 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                        {
                            cmplx da, db;
                            {
                                cmplx ca, cb; ca.r = t1.r + tw5r * t2.r + tw1r * t3.r + tw4r * t4.r + tw2r * t5.r + tw3r * t6.r; ca.i = t1.i + tw5r * t2.i + tw1r * t3.i + tw4r * t4.i + tw2r * t5.i + tw3r * t6.i; cb.i = +tw5i * t11.r - tw1i * t10.r + tw4i * t9.r - tw2i * t8.r + tw3i * t7.r; cb.r = -(+tw5i * t11.i - tw1i * t10.i + tw4i * t9.i - tw2i * t8.i + tw3i * t7.i);
                                {
                                    da.r = ca.r + cb.r; da.i = ca.i + cb.i; db.r = ca.r - cb.r; db.i = ca.i - cb.i;
                                }
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (5))].r = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * da.r - sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * da.i; ch[(i) + ido * ((k) + l1 * (5))].i = wa[(i) - 1 + (5 - 1) * (ido - 1)].r * da.i + sign * wa[(i) - 1 + (5 - 1) * (ido - 1)].i * da.r;
                            }
                            {
                                ch[(i) + ido * ((k) + l1 * (6))].r = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.r - sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.i; ch[(i) + ido * ((k) + l1 * (6))].i = wa[(i) - 1 + (6 - 1) * (ido - 1)].r * db.i + sign * wa[(i) - 1 + (6 - 1) * (ido - 1)].i * db.r;
                            }
                        }
                    }
                }
        }

        static void passg(int ido, int ip, int l1, Span<cmplx> cc, Span<cmplx> ch, Span<cmplx> wa, Span<cmplx> csarr, int sign)
        {
            int cdim = ip;
            int ipph = (ip + 1) / 2;
            int idl1 = ido * l1;

            if (NAIL_TEST) Console.WriteLine("start passg");
            if (NAIL_TEST) NailTestPrintComplexArray(cc.Slice(0, cdim * l1));
            Span<cmplx> wal = new cmplx[ip];
            wal[0] = new cmplx(1.0, 0.0);

            for (int i = 1; i < ip; ++i)
            {
                wal[i] = new cmplx(csarr[i].r, sign * csarr[i].i);
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 0; i < ido; ++i)
                {
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))];
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
            {
                for (int k = 0; k < l1; ++k)
                {
                    for (int i = 0; i < ido; ++i)
                    {
                        ch[(i) + ido * ((k) + l1 * (j))].r = cc[(i) + ido * ((j) + cdim * (k))].r + cc[(i) + ido * ((jc) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (j))].i = cc[(i) + ido * ((j) + cdim * (k))].i + cc[(i) + ido * ((jc) + cdim * (k))].i;
                        ch[(i) + ido * ((k) + l1 * (jc))].r = cc[(i) + ido * ((j) + cdim * (k))].r - cc[(i) + ido * ((jc) + cdim * (k))].r;
                        ch[(i) + ido * ((k) + l1 * (jc))].i = cc[(i) + ido * ((j) + cdim * (k))].i - cc[(i) + ido * ((jc) + cdim * (k))].i;
                    }
                }
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 0; i < ido; ++i)
                {
                    cmplx tmp = ch[(i) + ido * ((k) + l1 * (0))];
                    for (int j = 1; j < ipph; ++j)
                    {
                        tmp.r = tmp.r + ch[(i) + ido * ((k) + l1 * (j))].r;
                        tmp.i = tmp.i + ch[(i) + ido * ((k) + l1 * (j))].i;
                    }

                    cc[(i) + ido * ((k) + l1 * (0))] = tmp;
                }
            }

            for (int l = 1, lc = ip - 1; l < ipph; ++l, --lc)
            {
                // j=0
                for (int ik = 0; ik < idl1; ++ik)
                {
                    cc[(ik) + idl1 * (l)].r = ch[(ik) + idl1 * (0)].r + wal[l].r * ch[(ik) + idl1 * (1)].r + wal[2 * l].r * ch[(ik) + idl1 * (2)].r;
                    cc[(ik) + idl1 * (l)].i = ch[(ik) + idl1 * (0)].i + wal[l].r * ch[(ik) + idl1 * (1)].i + wal[2 * l].r * ch[(ik) + idl1 * (2)].i;
                    cc[(ik) + idl1 * (lc)].r = -wal[l].i * ch[(ik) + idl1 * (ip - 1)].i - wal[2 * l].i * ch[(ik) + idl1 * (ip - 2)].i;
                    cc[(ik) + idl1 * (lc)].i = wal[l].i * ch[(ik) + idl1 * (ip - 1)].r + wal[2 * l].i * ch[(ik) + idl1 * (ip - 2)].r;
                }

                int iwal = 2 * l;
                int j = 3, jc = ip - 3;
                for (; j < ipph - 1; j += 2, jc -= 2)
                {
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal = wal[iwal];
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal2 = wal[iwal];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)].r += ch[(ik) + idl1 * (j)].r * xwal.r + ch[(ik) + idl1 * (j + 1)].r * xwal2.r;
                        cc[(ik) + idl1 * (l)].i += ch[(ik) + idl1 * (j)].i * xwal.r + ch[(ik) + idl1 * (j + 1)].i * xwal2.r;
                        cc[(ik) + idl1 * (lc)].r -= ch[(ik) + idl1 * (jc)].i * xwal.i + ch[(ik) + idl1 * (jc - 1)].i * xwal2.i;
                        cc[(ik) + idl1 * (lc)].i += ch[(ik) + idl1 * (jc)].r * xwal.i + ch[(ik) + idl1 * (jc - 1)].r * xwal2.i;
                    }
                }

                for (; j < ipph; ++j, --jc)
                {
                    iwal += l; if (iwal > ip) iwal -= ip;
                    cmplx xwal = wal[iwal];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)].r += ch[(ik) + idl1 * (j)].r * xwal.r;
                        cc[(ik) + idl1 * (l)].i += ch[(ik) + idl1 * (j)].i * xwal.r;
                        cc[(ik) + idl1 * (lc)].r -= ch[(ik) + idl1 * (jc)].i * xwal.i;
                        cc[(ik) + idl1 * (lc)].i += ch[(ik) + idl1 * (jc)].r * xwal.i;
                    }
                }
            }
            //DEALLOC(wal);

            // shuffling and twiddling
            if (ido == 1)
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                {
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cmplx t1 = cc[(ik) + idl1 * (j)], t2 = cc[(ik) + idl1 * (jc)];
                        cc[(ik) + idl1 * (j)].r = t1.r + t2.r; cc[(ik) + idl1 * (j)].i = t1.i + t2.i;
                        cc[(ik) + idl1 * (jc)].r = t1.r - t2.r; cc[(ik) + idl1 * (jc)].i = t1.i - t2.i;
                    }
                }
            }
            else
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)
                {
                    for (int k = 0; k < l1; ++k)
                    {
                        cmplx t1 = cc[(0) + ido * ((k) + l1 * (j))], t2 = cc[(0) + ido * ((k) + l1 * (jc))];
                        cc[(0) + ido * ((k) + l1 * (j))].r = t1.r + t2.r;
                        cc[(0) + ido * ((k) + l1 * (j))].i = t1.i + t2.i;
                        cc[(0) + ido * ((k) + l1 * (jc))].r = t1.r - t2.r;
                        cc[(0) + ido * ((k) + l1 * (jc))].i = t1.i - t2.i;

                        for (int i = 1; i < ido; ++i)
                        {
                            cmplx x1, x2;
                            x1.r = cc[(i) + ido * ((k) + l1 * (j))].r + cc[(i) + ido * ((k) + l1 * (jc))].r;
                            x1.i = cc[(i) + ido * ((k) + l1 * (j))].i + cc[(i) + ido * ((k) + l1 * (jc))].i;
                            x2.r = cc[(i) + ido * ((k) + l1 * (j))].r - cc[(i) + ido * ((k) + l1 * (jc))].r;
                            x2.i = cc[(i) + ido * ((k) + l1 * (j))].i - cc[(i) + ido * ((k) + l1 * (jc))].i;

                            int idij = (j - 1) * (ido - 1) + i - 1;
                            cc[(i) + ido * ((k) + l1 * (j))].r = wa[idij].r * x1.r - sign * wa[idij].i * x1.i;
                            cc[(i) + ido * ((k) + l1 * (j))].i = wa[idij].r * x1.i + sign * wa[idij].i * x1.r;
                            idij = (jc - 1) * (ido - 1) + i - 1;
                            cc[(i) + ido * ((k) + l1 * (jc))].r = wa[idij].r * x2.r - sign * wa[idij].i * x2.i;
                            cc[(i) + ido * ((k) + l1 * (jc))].i = wa[idij].r * x2.i + sign * wa[idij].i * x2.r;
                        }
                    }
                }
            }
        }

        static void pass_all(cfftp_plan plan, Span<cmplx> c, double fct, int sign)
        {
            if (plan.length == 1)
            {
                return;
            }

            int len = plan.length;
            int l1 = 1, nf = plan.nfct;
            Span<cmplx> ch = new cmplx[len];
            Span<cmplx> p1 = c, p2 = ch;

            for (int k1 = 0; k1 < nf; k1++)
            {
                int ip = plan.fct[k1].fct;
                int l2 = ip * l1;
                int ido = len / l2;
                if (ip == 4)
                {
                    if (sign > 0)
                    {
                        pass4b(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                    else
                    {
                        pass4f(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                }
                else if (ip == 2)
                {
                    if (sign > 0)
                    {
                        pass2b(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                    else
                    {
                        pass2f(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                }
                else if (ip == 3)
                {
                    if (sign > 0)
                    {
                        pass3b(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                    else
                    {
                        pass3f(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                }
                else if (ip == 5)
                {
                    if (sign > 0)
                    {
                        pass5b(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                    else
                    {
                        pass5f(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw));
                    }
                }
                else if (ip == 7)
                {
                    pass7(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw), sign);
                }
                else if (ip == 11)
                {
                    pass11(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw), sign);
                }
                else
                {
                    passg(ido, ip, l1, p1, p2, plan.mem.AsSpan(plan.fct[k1].tw), plan.mem.AsSpan(plan.fct[k1].tws), sign);
                    {
                        Span<cmplx> tmp = p1;
                        p1 = p2;
                        p2 = tmp;
                    }
                }

                {
                    Span<cmplx> tmp = p1;
                    p1 = p2;
                    p2 = tmp;
                }
                l1 = l2;
            }

            if (p1 != c) // FIXME span address comparison?
            {
                if (fct != 1.0)
                    for (int i = 0; i < len; ++i)
                    {
                        c[i].r = ch[i].r * fct;
                        c[i].i = ch[i].i * fct;
                    }
                else
                {
                    p1.Slice(0, len).CopyTo(c);
                    //memcpy(c, p1, len * sizeof(cmplx));
                }
            }
            else
            {
                if (fct != 1.0)
                {
                    for (int i = 0; i < len; ++i)
                    {
                        c[i].r *= fct;
                        c[i].i *= fct;
                    }
                }
            }

            //DEALLOC(ch);
        }

        static void cfftp_forward(cfftp_plan plan, Span<double> c, double fct)
        {
            pass_all(plan, MemoryMarshal.Cast<double, cmplx>(c), fct, -1);
        }

        static void cfftp_backward(cfftp_plan plan, Span<double> c, double fct)
        {
            pass_all(plan, MemoryMarshal.Cast<double, cmplx>(c), fct, 1);
        }

        static bool cfftp_factorize(cfftp_plan plan)
        {
            int length = plan.length;
            int nfct = 0;

            while ((length % 4) == 0)
            {
                if (nfct >= NFCT)
                {
                    return false;
                }

                plan.fct[nfct++].fct = 4;
                length >>= 2;
            }
            if ((length % 2) == 0)
            {
                length >>= 1;
                // factor 2 should be at the front of the factor list
                if (nfct >= NFCT)
                {
                    return false;
                }

                plan.fct[nfct++].fct = 2;
                Swap(ref plan.fct[0].fct, ref plan.fct[nfct - 1].fct);
            }

            int maxl = (int)(Math.Sqrt((double)length)) + 1;
            for (int divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
            {
                if ((length % divisor) == 0)
                {
                    while ((length % divisor) == 0)
                    {
                        if (nfct >= NFCT)
                        {
                            return false;
                        }

                        plan.fct[nfct++].fct = divisor;
                        length /= divisor;
                    }

                    maxl = (int)(Math.Sqrt((double)length)) + 1;
                }
            }

            if (length > 1)
            {
                plan.fct[nfct++].fct = length;
            }

            plan.nfct = nfct;
            return true;
        }

        static int cfftp_twsize(cfftp_plan plan)
        {
            int twsize = 0, l1 = 1;
            for (int k = 0; k < plan.nfct; ++k)
            {
                int ip = plan.fct[k].fct, ido = plan.length / (l1 * ip);
                twsize += (ip - 1) * (ido - 1);
                if (ip > 11)
                {
                    twsize += ip;
                }

                l1 *= ip;
            }

            return twsize;
        }

        static void cfftp_comp_twiddle(cfftp_plan plan)
        {
            int length = plan.length;
            double[] twid = new double[2 * length];
            sincos_2pibyn(length, twid);
            int l1 = 1;
            int memofs = 0;
            for (int k = 0; k < plan.nfct; ++k)
            {
                int ip = plan.fct[k].fct, ido = length / (l1 * ip);
                plan.fct[k].tw = memofs;
                Span<cmplx> tw = plan.mem.AsSpan(memofs);
                memofs += (ip - 1) * (ido - 1);
                for (int j = 1; j < ip; ++j)
                {
                    for (int i = 1; i < ido; ++i)
                    {
                        tw[(j - 1) * (ido - 1) + i - 1].r = twid[2 * j * l1 * i];
                        tw[(j - 1) * (ido - 1) + i - 1].i = twid[2 * j * l1 * i + 1];
                    }
                }
                if (ip > 11)
                {
                    plan.fct[k].tws = memofs;
                    Span<cmplx> tws = plan.mem.AsSpan(memofs);
                    memofs += ip;
                    for (int j = 0; j < ip; ++j)
                    {
                        tws[j].r = twid[2 * j * l1 * ido];
                        tws[j].i = twid[2 * j * l1 * ido + 1];
                    }
                }
                l1 *= ip;
            }

            //DEALLOC(twid);
        }

        static cfftp_plan make_cfftp_plan(int length)
        {
            if (length == 0)
            {
                return null;
            }

            cfftp_plan plan = new cfftp_plan();
            plan.length = length;
            plan.nfct = 0;
#if !NET9_0_OR_GREATER
            if (plan.fct == null) // if it's not an inline struct array...
            {
                plan.fct = new cfftp_fctdata[NFCT];
            }
#endif

            for (int i = 0; i < NFCT; ++i)
            {
                plan.fct[i] = new cfftp_fctdata(0, 0, 0);
            }

            plan.mem = null;
            if (length == 1)
            {
                return plan;
            }

            if (!cfftp_factorize(plan))
            {
                return null;
            }

            int tws = cfftp_twsize(plan);
            plan.mem = new cmplx[tws];

            cfftp_comp_twiddle(plan);

            //if (NAIL_TEST) Console.WriteLine("cfftp plan mem");
            //if (NAIL_TEST) NailTestPrintComplexArray(plan.mem.AsSpan(0, tws));

            return plan;
        }

        static void destroy_cfftp_plan(cfftp_plan plan)
        {
        }

        public struct rfftp_fctdata
        {
            public int fct;
            public int tw; // index into plan.mem
            public int tws; // index into plan.mem

            public rfftp_fctdata(int fct, int tw, int tws)
            {
                this.fct = fct;
                this.tw = tw;
                this.tws = tws;
            }
        }

#if NET9_0_OR_GREATER
        [System.Runtime.CompilerServices.InlineArray(NFCT)]
        public struct rfftp_fctdata_array
        {
            public rfftp_fctdata data;
        }

        public class rfftp_plan
        {
            public int length;
            public int nfct;
            public double[] mem;
            public rfftp_fctdata_array fct;
        }
#else
        public class rfftp_plan
        {
            public int length;
            public int nfct;
            public double[] mem;
            public rfftp_fctdata[] fct; // [NFCT]

            public IEnumerable<int> Factors()
            {
                foreach (var data in fct)
                {
                    if (data.fct == 0)
                    {
                        yield break;
                    }

                    yield return data.fct;
                }
            }
        }
#endif


        static void radf2(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 2;

            if (NAIL_TEST) Console.WriteLine("start radf2");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + cc[(0) + ido * ((k) + l1 * (1))];
                ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] - cc[(0) + ido * ((k) + l1 * (1))];
            }

            if ((ido & 1) == 0)
            {
                for (int k = 0; k < l1; k++)
                {
                    ch[(0) + ido * ((1) + cdim * (k))] = -cc[(ido - 1) + ido * ((k) + l1 * (1))];
                    ch[(ido - 1) + ido * ((0) + cdim * (k))] = cc[(ido - 1) + ido * ((k) + l1 * (0))];
                }
            }

            if (ido <= 2)
            {
                return;
            }

            for (int k = 0; k < l1; k++)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double tr2, ti2;
                    tr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))];
                    ti2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))];
                    ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k) + l1 * (0))] + tr2;
                    ch[(ic - 1) + ido * ((1) + cdim * (k))] = cc[(i - 1) + ido * ((k) + l1 * (0))] - tr2;
                    ch[(i) + ido * ((0) + cdim * (k))] = ti2 + cc[(i) + ido * ((k) + l1 * (0))];
                    ch[(ic) + ido * ((1) + cdim * (k))] = ti2 - cc[(i) + ido * ((k) + l1 * (0))];
                }
            }
        }

        static void radf3(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 3;
            const double taur = -0.5, taui = 0.86602540378443864676;

            if (NAIL_TEST) Console.WriteLine("start radf3");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double cr2 = cc[(0) + ido * ((k) + l1 * (1))] + cc[(0) + ido * ((k) + l1 * (2))];
                ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + cr2;
                ch[(0) + ido * ((2) + cdim * (k))] = taui * (cc[(0) + ido * ((k) + l1 * (2))] - cc[(0) + ido * ((k) + l1 * (1))]);
                ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + taur * cr2;
            }

            if (ido == 1)
            {
                return;
            }

            for (int k = 0; k < l1; k++)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double di2, di3, dr2, dr3;
                    dr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))];
                    di2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))]; // d2=conj(WA0)*CC1
                    dr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))];
                    di3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))]; // d3=conj(WA1)*CC2
                    double cr2 = dr2 + dr3; // c add
                    double ci2 = di2 + di3;
                    ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k) + l1 * (0))] + cr2; // c add
                    ch[(i) + ido * ((0) + cdim * (k))] = cc[(i) + ido * ((k) + l1 * (0))] + ci2;
                    double tr2 = cc[(i - 1) + ido * ((k) + l1 * (0))] + taur * cr2; // c add
                    double ti2 = cc[(i) + ido * ((k) + l1 * (0))] + taur * ci2;
                    double tr3 = taui * (di2 - di3);  // t3 = taui*i*(d3-d2)?
                    double ti3 = taui * (dr3 - dr2);
                    ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr3;
                    ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr3; // PM(i) = t2+t3
                    ch[(i) + ido * ((2) + cdim * (k))] = ti3 + ti2;
                    ch[(ic) + ido * ((1) + cdim * (k))] = ti3 - ti2; // PM(ic) = conj(t2-t3)
                }
            }
        }

        static void radf4(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 4;
            const double hsqt2 = 0.70710678118654752440;

            if (NAIL_TEST) Console.WriteLine("start radf4");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double tr1, tr2;
                tr1 = cc[(0) + ido * ((k) + l1 * (3))] + cc[(0) + ido * ((k) + l1 * (1))];
                ch[(0) + ido * ((2) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (3))] - cc[(0) + ido * ((k) + l1 * (1))];
                tr2 = cc[(0) + ido * ((k) + l1 * (0))] + cc[(0) + ido * ((k) + l1 * (2))];
                ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] - cc[(0) + ido * ((k) + l1 * (2))];
                ch[(0) + ido * ((0) + cdim * (k))] = tr2 + tr1;
                ch[(ido - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
            }

            if ((ido & 1) == 0)
            {
                for (int k = 0; k < l1; k++)
                {
                    double ti1 = -hsqt2 * (cc[(ido - 1) + ido * ((k) + l1 * (1))] + cc[(ido - 1) + ido * ((k) + l1 * (3))]);
                    double tr1 = hsqt2 * (cc[(ido - 1) + ido * ((k) + l1 * (1))] - cc[(ido - 1) + ido * ((k) + l1 * (3))]);
                    ch[(ido - 1) + ido * ((0) + cdim * (k))] = cc[(ido - 1) + ido * ((k) + l1 * (0))] + tr1;
                    ch[(ido - 1) + ido * ((2) + cdim * (k))] = cc[(ido - 1) + ido * ((k) + l1 * (0))] - tr1;
                    ch[(0) + ido * ((3) + cdim * (k))] = ti1 + cc[(ido - 1) + ido * ((k) + l1 * (2))];
                    ch[(0) + ido * ((1) + cdim * (k))] = ti1 - cc[(ido - 1) + ido * ((k) + l1 * (2))];
                }
            }

            if (ido <= 2)
            {
                return;
            }

            for (int k = 0; k < l1; k++)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
                    cr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))];
                    ci2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))];
                    cr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))];
                    ci3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))];
                    cr4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (3))] + wa[(i - 1) + (2) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (3))];
                    ci4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (3))] - wa[(i - 1) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (3))];
                    tr1 = cr4 + cr2; tr4 = cr4 - cr2;
                    ti1 = ci2 + ci4; ti4 = ci2 - ci4;
                    tr2 = cc[(i - 1) + ido * ((k) + l1 * (0))] + cr3;
                    tr3 = cc[(i - 1) + ido * ((k) + l1 * (0))] - cr3;
                    ti2 = cc[(i) + ido * ((k) + l1 * (0))] + ci3;
                    ti3 = cc[(i) + ido * ((k) + l1 * (0))] - ci3;
                    ch[(i - 1) + ido * ((0) + cdim * (k))] = tr2 + tr1; ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
                    ch[(i) + ido * ((0) + cdim * (k))] = ti1 + ti2;
                    ch[(ic) + ido * ((3) + cdim * (k))] = ti1 - ti2;
                    ch[(i - 1) + ido * ((2) + cdim * (k))] = tr3 + ti4; ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr3 - ti4;
                    ch[(i) + ido * ((2) + cdim * (k))] = tr4 + ti3;
                    ch[(ic) + ido * ((1) + cdim * (k))] = tr4 - ti3;
                }
            }
        }

        static void radf5(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 5;
            const double 
                tr11 = 0.3090169943749474241,
                ti11 = 0.95105651629515357212,
                tr12 = -0.8090169943749474241,
                ti12 = 0.58778525229247312917;

            if (NAIL_TEST) Console.WriteLine("start radf5");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double cr2, cr3, ci4, ci5;
                cr2 = cc[(0) + ido * ((k) + l1 * (4))] + cc[(0) + ido * ((k) + l1 * (1))];
                ci5 = cc[(0) + ido * ((k) + l1 * (4))] - cc[(0) + ido * ((k) + l1 * (1))];
                cr3 = cc[(0) + ido * ((k) + l1 * (3))] + cc[(0) + ido * ((k) + l1 * (2))];
                ci4 = cc[(0) + ido * ((k) + l1 * (3))] - cc[(0) + ido * ((k) + l1 * (2))];
                ch[(0) + ido * ((0) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + cr2 + cr3;
                ch[(ido - 1) + ido * ((1) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + tr11 * cr2 + tr12 * cr3;
                ch[(0) + ido * ((2) + cdim * (k))] = ti11 * ci5 + ti12 * ci4;
                ch[(ido - 1) + ido * ((3) + cdim * (k))] = cc[(0) + ido * ((k) + l1 * (0))] + tr12 * cr2 + tr11 * cr3;
                ch[(0) + ido * ((4) + cdim * (k))] = ti12 * ci5 - ti11 * ci4;
            }

            if (ido == 1)
            {
                return;
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    double ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3,
                        dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
                    int ic = ido - i;
                    dr2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))] + wa[(i - 1) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))];
                    di2 = wa[(i - 2) + (0) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (1))] - wa[(i - 1) + (0) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (1))];
                    dr3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))] + wa[(i - 1) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))];
                    di3 = wa[(i - 2) + (1) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (2))] - wa[(i - 1) + (1) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (2))];
                    dr4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (3))] + wa[(i - 1) + (2) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (3))];
                    di4 = wa[(i - 2) + (2) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (3))] - wa[(i - 1) + (2) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (3))];
                    dr5 = wa[(i - 2) + (3) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (4))] + wa[(i - 1) + (3) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (4))];
                    di5 = wa[(i - 2) + (3) * (ido - 1)] * cc[(i) + ido * ((k) + l1 * (4))] - wa[(i - 1) + (3) * (ido - 1)] * cc[(i - 1) + ido * ((k) + l1 * (4))];
                    cr2 = dr5 + dr2; ci5 = dr5 - dr2;
                    ci2 = di2 + di5; cr5 = di2 - di5;
                    cr3 = dr4 + dr3; ci4 = dr4 - dr3;
                    ci3 = di3 + di4; cr4 = di3 - di4;
                    ch[(i - 1) + ido * ((0) + cdim * (k))] = cc[(i - 1) + ido * ((k) + l1 * (0))] + cr2 + cr3;
                    ch[(i) + ido * ((0) + cdim * (k))] = cc[(i) + ido * ((k) + l1 * (0))] + ci2 + ci3;
                    tr2 = cc[(i - 1) + ido * ((k) + l1 * (0))] + tr11 * cr2 + tr12 * cr3;
                    ti2 = cc[(i) + ido * ((k) + l1 * (0))] + tr11 * ci2 + tr12 * ci3;
                    tr3 = cc[(i - 1) + ido * ((k) + l1 * (0))] + tr12 * cr2 + tr11 * cr3;
                    ti3 = cc[(i) + ido * ((k) + l1 * (0))] + tr12 * ci2 + tr11 * ci3;
                    tr5 = cr5 * ti11 + cr4 * ti12;
                    tr4 = cr5 * ti12 - cr4 * ti11;
                    ti5 = ci5 * ti11 + ci4 * ti12;
                    ti4 = ci5 * ti12 - ci4 * ti11;
                    ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr5;
                    ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr5;
                    ch[(i) + ido * ((2) + cdim * (k))] = ti5 + ti2;
                    ch[(ic) + ido * ((1) + cdim * (k))] = ti5 - ti2;
                    ch[(i - 1) + ido * ((4) + cdim * (k))] = tr3 + tr4;
                    ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr3 - tr4;
                    ch[(i) + ido * ((4) + cdim * (k))] = ti4 + ti3;
                    ch[(ic) + ido * ((3) + cdim * (k))] = ti4 - ti3;
                }
            }
        }

        static void radfg(int ido, int ip, int l1, Span<double> cc, Span<double> ch, Span<double> wa, Span<double> csarr)
        {
            int cdim = ip;
            int ipph = (ip + 1) / 2;
            int idl1 = ido * l1;

            if (NAIL_TEST) Console.WriteLine("start radfg");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            if (ido > 1)
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)              // 114
                {
                    int is1 = (j - 1) * (ido - 1),
                        is2 = (jc - 1) * (ido - 1);
                    for (int k = 0; k < l1; ++k)                            // 113
                    {
                        int idij = is1;
                        int idij2 = is2;

                        for (int i = 1; i <= ido - 2; i += 2)                      // 112
                        {
                            double t1 = cc[(i) + ido * ((k) + l1 * (j))], t2 = cc[(i + 1) + ido * ((k) + l1 * (j))],
                                t3 = cc[(i) + ido * ((k) + l1 * (jc))], t4 = cc[(i + 1) + ido * ((k) + l1 * (jc))];
                            double x1 = wa[idij] * t1 + wa[idij + 1] * t2,
                                x2 = wa[idij] * t2 - wa[idij + 1] * t1,
                                x3 = wa[idij2] * t3 + wa[idij2 + 1] * t4,
                                x4 = wa[idij2] * t4 - wa[idij2 + 1] * t3;
                            cc[(i) + ido * ((k) + l1 * (j))] = x1 + x3;
                            cc[(i) + ido * ((k) + l1 * (jc))] = x2 - x4;
                            cc[(i + 1) + ido * ((k) + l1 * (j))] = x2 + x4;
                            cc[(i + 1) + ido * ((k) + l1 * (jc))] = x3 - x1;
                            idij += 2;
                            idij2 += 2;
                        }
                    }
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 123
            {
                for (int k = 0; k < l1; ++k)                              // 122
                {
                    double t1 = cc[(0) + ido * ((k) + l1 * (j))], t2 = cc[(0) + ido * ((k) + l1 * (jc))];
                    cc[(0) + ido * ((k) + l1 * (j))] = t1 + t2;
                    cc[(0) + ido * ((k) + l1 * (jc))] = t2 - t1;
                }
            }

            //everything in C
            //memset(ch,0,ip*l1*ido*sizeof(double));

            for (int l = 1, lc = ip - 1; l < ipph; ++l, --lc)                 // 127
            {
                for (int ik = 0; ik < idl1; ++ik)                         // 124
                {
                    ch[(ik) + idl1 * (l)] = cc[(ik) + idl1 * (0)] + csarr[2 * l] * cc[(ik) + idl1 * (1)] + csarr[4 * l] * cc[(ik) + idl1 * (2)];
                    ch[(ik) + idl1 * (lc)] = csarr[2 * l + 1] * cc[(ik) + idl1 * (ip - 1)] + csarr[4 * l + 1] * cc[(ik) + idl1 * (ip - 2)];
                }

                int iang = 2 * l;
                int j = 3, jc = ip - 3;
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
                    for (int ik = 0; ik < idl1; ++ik)                       // 125
                    {
                        ch[(ik) + idl1 * (l)] += ar1 * cc[(ik) + idl1 * (j)] + ar2 * cc[(ik) + idl1 * (j + 1)]
                            + ar3 * cc[(ik) + idl1 * (j + 2)] + ar4 * cc[(ik) + idl1 * (j + 3)];
                        ch[(ik) + idl1 * (lc)] += ai1 * cc[(ik) + idl1 * (jc)] + ai2 * cc[(ik) + idl1 * (jc - 1)]
                            + ai3 * cc[(ik) + idl1 * (jc - 2)] + ai4 * cc[(ik) + idl1 * (jc - 3)];
                    }
                }

                for (; j < ipph - 1; j += 2, jc -= 2)              // 126
                {
                    iang += l; if (iang >= ip) iang -= ip;
                    double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
                    iang += l; if (iang >= ip) iang -= ip;
                    double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
                    for (int ik = 0; ik < idl1; ++ik)                       // 125
                    {
                        ch[(ik) + idl1 * (l)] += ar1 * cc[(ik) + idl1 * (j)] + ar2 * cc[(ik) + idl1 * (j + 1)];
                        ch[(ik) + idl1 * (lc)] += ai1 * cc[(ik) + idl1 * (jc)] + ai2 * cc[(ik) + idl1 * (jc - 1)];
                    }
                }

                for (; j < ipph; ++j, --jc)              // 126
                {
                    iang += l; if (iang >= ip) iang -= ip;
                    double ar = csarr[2 * iang], ai = csarr[2 * iang + 1];
                    for (int ik = 0; ik < idl1; ++ik)                       // 125
                    {
                        ch[(ik) + idl1 * (l)] += ar * cc[(ik) + idl1 * (j)];
                        ch[(ik) + idl1 * (lc)] += ai * cc[(ik) + idl1 * (jc)];
                    }
                }
            }
            for (int ik = 0; ik < idl1; ++ik)                         // 101
            {
                ch[(ik) + idl1 * (0)] = cc[(ik) + idl1 * (0)];
            }

            for (int j = 1; j < ipph; ++j)                              // 129
            {
                for (int ik = 0; ik < idl1; ++ik)                         // 128
                {
                    ch[(ik) + idl1 * (0)] += cc[(ik) + idl1 * (j)];
                }
            }

            // everything in CH at this point!
            //memset(cc,0,ip*l1*ido*sizeof(double));

            for (int k = 0; k < l1; ++k)                                // 131
            {
                for (int i = 0; i < ido; ++i)                             // 130
                {
                    cc[(i) + ido * ((0) + cdim * (k))] = ch[(i) + ido * ((k) + l1 * (0))];
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 137
            {
                int j2 = 2 * j - 1;
                for (int k = 0; k < l1; ++k)                              // 136
                {
                    cc[(ido - 1) + ido * ((j2) + cdim * (k))] = ch[(0) + ido * ((k) + l1 * (j))];
                    cc[(0) + ido * ((j2 + 1) + cdim * (k))] = ch[(0) + ido * ((k) + l1 * (jc))];
                }
            }

            if (ido == 1)
            {
                return;
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)                // 140
            {
                int j2 = 2 * j - 1;
                for (int k = 0; k < l1; ++k)                               // 139
                {
                    for (int i = 1, ic = ido - i - 2; i <= ido - 2; i += 2, ic -= 2)      // 138
                    {
                        cc[(i) + ido * ((j2 + 1) + cdim * (k))] = ch[(i) + ido * ((k) + l1 * (j))] + ch[(i) + ido * ((k) + l1 * (jc))];
                        cc[(ic) + ido * ((j2) + cdim * (k))] = ch[(i) + ido * ((k) + l1 * (j))] - ch[(i) + ido * ((k) + l1 * (jc))];
                        cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] = ch[(i + 1) + ido * ((k) + l1 * (j))] + ch[(i + 1) + ido * ((k) + l1 * (jc))];
                        cc[(ic + 1) + ido * ((j2) + cdim * (k))] = ch[(i + 1) + ido * ((k) + l1 * (jc))] - ch[(i + 1) + ido * ((k) + l1 * (j))];
                    }
                }
            }
        }

        static void radb2(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 2;

            if (NAIL_TEST) Console.WriteLine("start radb2");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                ch[(0) + ido * ((k) + l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((1) + cdim * (k))];
                ch[(0) + ido * ((k) + l1 * (1))] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((1) + cdim * (k))];
            }
            if ((ido & 1) == 0)
            {
                for (int k = 0; k < l1; k++)
                {
                    ch[(ido - 1) + ido * ((k) + l1 * (0))] = 2.0 * cc[(ido - 1) + ido * ((0) + cdim * (k))];
                    ch[(ido - 1) + ido * ((k) + l1 * (1))] = -2.0 * cc[(0) + ido * ((1) + cdim * (k))];
                }
            }

            if (ido <= 2)
            {
                return;
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double ti2, tr2;
                    ch[(i - 1) + ido * ((k) + l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))];
                    tr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))];
                    ti2 = cc[(i) + ido * ((0) + cdim * (k))] + cc[(ic) + ido * ((1) + cdim * (k))];
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))] - cc[(ic) + ido * ((1) + cdim * (k))];
                    ch[(i) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * ti2 + wa[(i - 1) + (0) * (ido - 1)] * tr2;
                    ch[(i - 1) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * tr2 - wa[(i - 1) + (0) * (ido - 1)] * ti2;
                }
            }
        }

        static void radb3(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 3;
            const double taur = -0.5, taui = 0.86602540378443864676;

            if (NAIL_TEST) Console.WriteLine("start radb3");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double tr2 = 2.0 * cc[(ido - 1) + ido * ((1) + cdim * (k))];
                double cr2 = cc[(0) + ido * ((0) + cdim * (k))] + taur * tr2;
                ch[(0) + ido * ((k) + l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + tr2;
                double ci3 = 2.0 * taui * cc[(0) + ido * ((2) + cdim * (k))];
                { ch[(0) + ido * ((k) + l1 * (2))] = cr2 + ci3; ch[(0) + ido * ((k) + l1 * (1))] = cr2 - ci3; };
            }

            if (ido == 1)
            {
                return;
            }

            for (int k = 0; k < l1; k++)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double tr2 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; // t2=CC(I) + conj(CC(ic))
                    double ti2 = cc[(i) + ido * ((2) + cdim * (k))] - cc[(ic) + ido * ((1) + cdim * (k))];
                    double cr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + taur * tr2;     // c2=CC +taur*t2
                    double ci2 = cc[(i) + ido * ((0) + cdim * (k))] + taur * ti2;
                    ch[(i - 1) + ido * ((k) + l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2;         // CH=CC+t2
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))] + ti2;
                    double cr3 = taui * (cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))]);// c3=taui*(CC(i)-conj(CC(ic)))
                    double ci3 = taui * (cc[(i) + ido * ((2) + cdim * (k))] + cc[(ic) + ido * ((1) + cdim * (k))]);
                    double dr3 = cr2 + ci3;
                    double dr2 = cr2 - ci3; // d2= (cr2-ci3, ci2+cr3) = c2+i*c3
                    double di2 = ci2 + cr3;
                    double di3 = ci2 - cr3; // d3= (cr2+ci3, ci2-cr3) = c2-i*c3
                    ch[(i) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * di2 + wa[(i - 1) + (0) * (ido - 1)] * dr2;
                    ch[(i - 1) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * dr2 - wa[(i - 1) + (0) * (ido - 1)] * di2; // ch = WA*d2
                    ch[(i) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * di3 + wa[(i - 1) + (1) * (ido - 1)] * dr3;
                    ch[(i - 1) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * dr3 - wa[(i - 1) + (1) * (ido - 1)] * di3;
                }
            }
        }

        static void radb4(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 4;
            const double sqrt2 = 1.41421356237309504880;

            if (NAIL_TEST) Console.WriteLine("start radb4");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double tr1, tr2;
                tr2 = cc[(0) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((3) + cdim * (k))];
                tr1 = cc[(0) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((3) + cdim * (k))];
                double tr3 = 2.0 * cc[(ido - 1) + ido * ((1) + cdim * (k))];
                double tr4 = 2.0 * cc[(0) + ido * ((2) + cdim * (k))];
                ch[(0) + ido * ((k) + l1 * (0))] = tr2 + tr3;
                ch[(0) + ido * ((k) + l1 * (2))] = tr2 - tr3;
                ch[(0) + ido * ((k) + l1 * (3))] = tr1 + tr4;
                ch[(0) + ido * ((k) + l1 * (1))] = tr1 - tr4;
            }
            if ((ido & 1) == 0)
            {
                for (int k = 0; k < l1; k++)
                {
                    double tr1, tr2, ti1, ti2;
                    ti1 = cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((1) + cdim * (k))];
                    ti2 = cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))];
                    tr2 = cc[(ido - 1) + ido * ((0) + cdim * (k))] + cc[(ido - 1) + ido * ((2) + cdim * (k))];
                    tr1 = cc[(ido - 1) + ido * ((0) + cdim * (k))] - cc[(ido - 1) + ido * ((2) + cdim * (k))];
                    ch[(ido - 1) + ido * ((k) + l1 * (0))] = tr2 + tr2;
                    ch[(ido - 1) + ido * ((k) + l1 * (1))] = sqrt2 * (tr1 - ti1);
                    ch[(ido - 1) + ido * ((k) + l1 * (2))] = ti2 + ti2;
                    ch[(ido - 1) + ido * ((k) + l1 * (3))] = -sqrt2 * (tr1 + ti1);
                }
            }

            if (ido <= 2)
            {
                return;
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
                    int ic = ido - i;
                    tr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + cc[(ic - 1) + ido * ((3) + cdim * (k))];
                    tr1 = cc[(i - 1) + ido * ((0) + cdim * (k))] - cc[(ic - 1) + ido * ((3) + cdim * (k))];
                    ti1 = cc[(i) + ido * ((0) + cdim * (k))] + cc[(ic) + ido * ((3) + cdim * (k))];
                    ti2 = cc[(i) + ido * ((0) + cdim * (k))] - cc[(ic) + ido * ((3) + cdim * (k))];
                    tr4 = cc[(i) + ido * ((2) + cdim * (k))] + cc[(ic) + ido * ((1) + cdim * (k))];
                    ti3 = cc[(i) + ido * ((2) + cdim * (k))] - cc[(ic) + ido * ((1) + cdim * (k))];
                    tr3 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))];
                    ti4 = cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))];
                    ch[(i - 1) + ido * ((k) + l1 * (0))] = tr2 + tr3; cr3 = tr2 - tr3;
                    ch[(i) + ido * ((k) + l1 * (0))] = ti2 + ti3; ci3 = ti2 - ti3;
                    cr4 = tr1 + tr4;
                    cr2 = tr1 - tr4;
                    ci2 = ti1 + ti4;
                    ci4 = ti1 - ti4;
                    ch[(i) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * ci2 + wa[(i - 1) + (0) * (ido - 1)] * cr2;
                    ch[(i - 1) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * cr2 - wa[(i - 1) + (0) * (ido - 1)] * ci2;
                    ch[(i) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * ci3 + wa[(i - 1) + (1) * (ido - 1)] * cr3;
                    ch[(i - 1) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * cr3 - wa[(i - 1) + (1) * (ido - 1)] * ci3;
                    ch[(i) + ido * ((k) + l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * ci4 + wa[(i - 1) + (2) * (ido - 1)] * cr4;
                    ch[(i - 1) + ido * ((k) + l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * cr4 - wa[(i - 1) + (2) * (ido - 1)] * ci4;
                }
            }
        }

        static void radb5(int ido, int l1, Span<double> cc, Span<double> ch, Span<double> wa)
        {
            const int cdim = 5;
            const double tr11 = 0.3090169943749474241, ti11 = 0.95105651629515357212,
                tr12 = -0.8090169943749474241, ti12 = 0.58778525229247312917;

            if (NAIL_TEST) Console.WriteLine("start radb5");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; k++)
            {
                double ti5 = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                double ti4 = cc[(0) + ido * ((4) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                double tr2 = cc[(ido - 1) + ido * ((1) + cdim * (k))] + cc[(ido - 1) + ido * ((1) + cdim * (k))];
                double tr3 = cc[(ido - 1) + ido * ((3) + cdim * (k))] + cc[(ido - 1) + ido * ((3) + cdim * (k))];
                ch[(0) + ido * ((k) + l1 * (0))] = cc[(0) + ido * ((0) + cdim * (k))] + tr2 + tr3;
                double cr2 = cc[(0) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
                double cr3 = cc[(0) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
                double ci4, ci5;
                ci5 = ti5 * ti11 + ti4 * ti12;
                ci4 = ti5 * ti12 - ti4 * ti11;
                ch[(0) + ido * ((k) + l1 * (4))] = cr2 + ci5;
                ch[(0) + ido * ((k) + l1 * (1))] = cr2 - ci5;
                ch[(0) + ido * ((k) + l1 * (3))] = cr3 + ci4;
                ch[(0) + ido * ((k) + l1 * (2))] = cr3 - ci4;
            }

            if (ido == 1)
            {
                return;
            }

            for (int k = 0; k < l1; ++k)
            {
                for (int i = 2; i < ido; i += 2)
                {
                    int ic = ido - i;
                    double tr2, tr3, tr4, tr5, ti2, ti3, ti4, ti5;
                    tr2 = cc[(i - 1) + ido * ((2) + cdim * (k))] + cc[(ic - 1) + ido * ((1) + cdim * (k))]; tr5 = cc[(i - 1) + ido * ((2) + cdim * (k))] - cc[(ic - 1) + ido * ((1) + cdim * (k))];
                    ti5 = cc[(i) + ido * ((2) + cdim * (k))] + cc[(ic) + ido * ((1) + cdim * (k))]; ti2 = cc[(i) + ido * ((2) + cdim * (k))] - cc[(ic) + ido * ((1) + cdim * (k))];
                    tr3 = cc[(i - 1) + ido * ((4) + cdim * (k))] + cc[(ic - 1) + ido * ((3) + cdim * (k))]; tr4 = cc[(i - 1) + ido * ((4) + cdim * (k))] - cc[(ic - 1) + ido * ((3) + cdim * (k))];
                    ti4 = cc[(i) + ido * ((4) + cdim * (k))] + cc[(ic) + ido * ((3) + cdim * (k))]; ti3 = cc[(i) + ido * ((4) + cdim * (k))] - cc[(ic) + ido * ((3) + cdim * (k))];
                    ch[(i - 1) + ido * ((k) + l1 * (0))] = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2 + tr3;
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))] + ti2 + ti3;
                    double cr2 = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
                    double ci2 = cc[(i) + ido * ((0) + cdim * (k))] + tr11 * ti2 + tr12 * ti3;
                    double cr3 = cc[(i - 1) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
                    double ci3 = cc[(i) + ido * ((0) + cdim * (k))] + tr12 * ti2 + tr11 * ti3;
                    double ci4, ci5, cr5, cr4;
                    cr5 = tr5 * ti11 + tr4 * ti12;
                    cr4 = tr5 * ti12 - tr4 * ti11;
                    ci5 = ti5 * ti11 + ti4 * ti12;
                    ci4 = ti5 * ti12 - ti4 * ti11;
                    double dr2, dr3, dr4, dr5, di2, di3, di4, di5;
                    dr4 = cr3 + ci4; dr3 = cr3 - ci4;
                    di3 = ci3 + cr4; di4 = ci3 - cr4;
                    dr5 = cr2 + ci5; dr2 = cr2 - ci5;
                    di2 = ci2 + cr5; di5 = ci2 - cr5;
                    ch[(i) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * di2 + wa[(i - 1) + (0) * (ido - 1)] * dr2;
                    ch[(i - 1) + ido * ((k) + l1 * (1))] = wa[(i - 2) + (0) * (ido - 1)] * dr2 - wa[(i - 1) + (0) * (ido - 1)] * di2;
                    ch[(i) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * di3 + wa[(i - 1) + (1) * (ido - 1)] * dr3;
                    ch[(i - 1) + ido * ((k) + l1 * (2))] = wa[(i - 2) + (1) * (ido - 1)] * dr3 - wa[(i - 1) + (1) * (ido - 1)] * di3;
                    ch[(i) + ido * ((k) + l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * di4 + wa[(i - 1) + (2) * (ido - 1)] * dr4;
                    ch[(i - 1) + ido * ((k) + l1 * (3))] = wa[(i - 2) + (2) * (ido - 1)] * dr4 - wa[(i - 1) + (2) * (ido - 1)] * di4;
                    ch[(i) + ido * ((k) + l1 * (4))] = wa[(i - 2) + (3) * (ido - 1)] * di5 + wa[(i - 1) + (3) * (ido - 1)] * dr5;
                    ch[(i - 1) + ido * ((k) + l1 * (4))] = wa[(i - 2) + (3) * (ido - 1)] * dr5 - wa[(i - 1) + (3) * (ido - 1)] * di5;
                }
            }
        }

        static void radbg(int ido, int ip, int l1, Span<double> cc, Span<double> ch, Span<double> wa, Span<double> csarr)
        {
            int cdim = ip;
            int ipph = (ip + 1) / 2;
            int idl1 = ido * l1;

            if (NAIL_TEST) Console.WriteLine("start radbg");
            if (NAIL_TEST) NailTestPrintDoubleArray(cc.Slice(0, cdim * l1));
            for (int k = 0; k < l1; ++k)        // 102
            {
                for (int i = 0; i < ido; ++i)     // 101
                {
                    ch[(i) + ido * ((k) + l1 * (0))] = cc[(i) + ido * ((0) + cdim * (k))];
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 108
            {
                int j2 = 2 * j - 1;
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (j))] = 2 * cc[(ido - 1) + ido * ((j2) + cdim * (k))];
                    ch[(0) + ido * ((k) + l1 * (jc))] = 2 * cc[(0) + ido * ((j2 + 1) + cdim * (k))];
                }
            }

            if (ido != 1)
            {
                for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 111
                {
                    int j2 = 2 * j - 1;
                    for (int k = 0; k < l1; ++k)
                    {
                        for (int i = 1, ic = ido - i - 2; i <= ido - 2; i += 2, ic -= 2)      // 109
                        {
                            ch[(i) + ido * ((k) + l1 * (j))] = cc[(i) + ido * ((j2 + 1) + cdim * (k))] + cc[(ic) + ido * ((j2) + cdim * (k))];
                            ch[(i) + ido * ((k) + l1 * (jc))] = cc[(i) + ido * ((j2 + 1) + cdim * (k))] - cc[(ic) + ido * ((j2) + cdim * (k))];
                            ch[(i + 1) + ido * ((k) + l1 * (j))] = cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] - cc[(ic + 1) + ido * ((j2) + cdim * (k))];
                            ch[(i + 1) + ido * ((k) + l1 * (jc))] = cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] + cc[(ic + 1) + ido * ((j2) + cdim * (k))];
                        }
                    }
                }
            }
            for (int l = 1, lc = ip - 1; l < ipph; ++l, --lc)
            {
                for (int ik = 0; ik < idl1; ++ik)
                {
                    cc[(ik) + idl1 * (l)] = ch[(ik) + idl1 * (0)] + csarr[2 * l] * ch[(ik) + idl1 * (1)] + csarr[4 * l] * ch[(ik) + idl1 * (2)];
                    cc[(ik) + idl1 * (lc)] = csarr[2 * l + 1] * ch[(ik) + idl1 * (ip - 1)] + csarr[4 * l + 1] * ch[(ik) + idl1 * (ip - 2)];
                }

                int iang = 2 * l;
                int j = 3, jc = ip - 3;
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
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)] += ar1 * ch[(ik) + idl1 * (j)] + ar2 * ch[(ik) + idl1 * (j + 1)]
                            + ar3 * ch[(ik) + idl1 * (j + 2)] + ar4 * ch[(ik) + idl1 * (j + 3)];
                        cc[(ik) + idl1 * (lc)] += ai1 * ch[(ik) + idl1 * (jc)] + ai2 * ch[(ik) + idl1 * (jc - 1)]
                            + ai3 * ch[(ik) + idl1 * (jc - 2)] + ai4 * ch[(ik) + idl1 * (jc - 3)];
                    }
                }
                for (; j < ipph - 1; j += 2, jc -= 2)
                {
                    iang += l; if (iang > ip) iang -= ip;
                    double ar1 = csarr[2 * iang], ai1 = csarr[2 * iang + 1];
                    iang += l; if (iang > ip) iang -= ip;
                    double ar2 = csarr[2 * iang], ai2 = csarr[2 * iang + 1];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)] += ar1 * ch[(ik) + idl1 * (j)] + ar2 * ch[(ik) + idl1 * (j + 1)];
                        cc[(ik) + idl1 * (lc)] += ai1 * ch[(ik) + idl1 * (jc)] + ai2 * ch[(ik) + idl1 * (jc - 1)];
                    }
                }
                for (; j < ipph; ++j, --jc)
                {
                    iang += l; if (iang > ip) iang -= ip;
                    double war = csarr[2 * iang], wai = csarr[2 * iang + 1];
                    for (int ik = 0; ik < idl1; ++ik)
                    {
                        cc[(ik) + idl1 * (l)] += war * ch[(ik) + idl1 * (j)];
                        cc[(ik) + idl1 * (lc)] += wai * ch[(ik) + idl1 * (jc)];
                    }
                }
            }

            for (int j = 1; j < ipph; ++j)
            {
                for (int ik = 0; ik < idl1; ++ik)
                {
                    ch[(ik) + idl1 * (0)] += ch[(ik) + idl1 * (j)];
                }
            }

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)   // 124
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (j))] = cc[(0) + ido * ((k) + l1 * (j))] - cc[(0) + ido * ((k) + l1 * (jc))];
                    ch[(0) + ido * ((k) + l1 * (jc))] = cc[(0) + ido * ((k) + l1 * (j))] + cc[(0) + ido * ((k) + l1 * (jc))];
                }
            }

            if (ido == 1) return;

            for (int j = 1, jc = ip - 1; j < ipph; ++j, --jc)  // 127
            {
                for (int k = 0; k < l1; ++k)
                {
                    for (int i = 1; i <= ido - 2; i += 2)
                    {
                        ch[(i) + ido * ((k) + l1 * (j))] = cc[(i) + ido * ((k) + l1 * (j))] - cc[(i + 1) + ido * ((k) + l1 * (jc))];
                        ch[(i) + ido * ((k) + l1 * (jc))] = cc[(i) + ido * ((k) + l1 * (j))] + cc[(i + 1) + ido * ((k) + l1 * (jc))];
                        ch[(i + 1) + ido * ((k) + l1 * (j))] = cc[(i + 1) + ido * ((k) + l1 * (j))] + cc[(i) + ido * ((k) + l1 * (jc))];
                        ch[(i + 1) + ido * ((k) + l1 * (jc))] = cc[(i + 1) + ido * ((k) + l1 * (j))] - cc[(i) + ido * ((k) + l1 * (jc))];
                    }
                }
            }

            // All in CH

            for (int j = 1; j < ip; ++j)
            {
                int is1 = (j - 1) * (ido - 1);
                for (int k = 0; k < l1; ++k)
                {
                    int idij = is1;
                    for (int i = 1; i <= ido - 2; i += 2)
                    {
                        double t1 = ch[(i) + ido * ((k) + l1 * (j))], t2 = ch[(i + 1) + ido * ((k) + l1 * (j))];
                        ch[(i) + ido * ((k) + l1 * (j))] = wa[idij] * t1 - wa[idij + 1] * t2;
                        ch[(i + 1) + ido * ((k) + l1 * (j))] = wa[idij] * t2 + wa[idij + 1] * t1;
                        idij += 2;
                    }
                }
            }
        }

        static void copy_and_norm(Span<double> c, Span<double> p1, int n, double fct)
        {
            if (p1 != c)
            {
                if (fct != 1.0)
                {
                    for (int i = 0; i < n; ++i)
                    {
                        c[i] = fct * p1[i];
                    }
                }
                else
                {
                    p1.Slice(0, n).CopyTo(c);
                    //memcpy(c, p1, n * sizeof(double));
                }
            }
            else
            {
                if (fct != 1.0)
                {
                    for (int i = 0; i < n; ++i)
                    {
                        c[i] *= fct;
                    }
                }
            }
        }

        static void rfftp_forward(rfftp_plan plan, Span<double> c, double fct)
        {
            if (plan.length == 1)
            {
                return;
            }

            int n = plan.length;
            int l1 = n, nf = plan.nfct;
            double[] ch = new double[n];
            Span<double> p1 = c, p2 = ch;

            for (int k1 = 0; k1 < nf; ++k1)
            {
                int k = nf - k1 - 1;
                int ip = plan.fct[k].fct;
                int ido = n / l1;
                l1 /= ip;

                if (ip == 4)
                {
                    radf4(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 2)
                {
                    radf2(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 3)
                {
                    radf3(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 5)
                {
                    radf5(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else
                {
                    radfg(ido, ip, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw), plan.mem.AsSpan(plan.fct[k].tws));
                    {
                        Span<double> tmp = p1;
                        p1 = p2;
                        p2 = tmp;
                    }
                }

                {
                    Span<double> tmp = p1;
                    p1 = p2;
                    p2 = tmp;
                }
            }

            copy_and_norm(c, p1, n, fct);
            //DEALLOC(ch);
        }

        static void rfftp_backward(rfftp_plan plan, Span<double> c, double fct)
        {
            if (plan.length == 1)
            {
                return;
            }

            int n = plan.length;
            int l1 = 1, nf = plan.nfct;
            double[] ch = new double[n];
            Span<double> p1 = c, p2 = ch;

            for (int k = 0; k < nf; k++)
            {
                int ip = plan.fct[k].fct,
                    ido = n / (ip * l1);

                if (ip == 4)
                {
                    radb4(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 2)
                {
                    radb2(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 3)
                {
                    radb3(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else if (ip == 5)
                {
                    radb5(ido, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw));
                }
                else
                {
                    radbg(ido, ip, l1, p1, p2, plan.mem.AsSpan(plan.fct[k].tw), plan.mem.AsSpan(plan.fct[k].tws));
                }

                {
                    Span<double> tmp = p1;
                    p1 = p2;
                    p2 = tmp;
                }
                l1 *= ip;
            }

            copy_and_norm(c, p1, n, fct);
            //DEALLOC(ch);
        }

        static bool rfftp_factorize(rfftp_plan plan)
        {
            int length = plan.length;
            int nfct = 0;
            while ((length % 4) == 0)
            {
                if (nfct >= NFCT)
                {
                    return false;
                }

                plan.fct[nfct++].fct = 4; length >>= 2;
            }

            if ((length % 2) == 0)
            {
                length >>= 1;
                // factor 2 should be at the front of the factor list
                if (nfct >= NFCT)
                {
                    return false;
                }

                plan.fct[nfct++].fct = 2;
                Swap(ref plan.fct[0].fct, ref plan.fct[nfct - 1].fct);
            }

            int maxl = (int)(Math.Sqrt((double)length)) + 1;
            for (int divisor = 3; (length > 1) && (divisor < maxl); divisor += 2)
            {
                if ((length % divisor) == 0)
                {
                    while ((length % divisor) == 0)
                    {
                        if (nfct >= NFCT)
                        {
                            return false;
                        }

                        plan.fct[nfct++].fct = divisor;
                        length /= divisor;
                    }
                    maxl = (int)(Math.Sqrt((double)length)) + 1;
                }
            }
            if (length > 1)
            {
                plan.fct[nfct++].fct = length;
            }

            plan.nfct = nfct;
            return true;
        }

        static int rfftp_twsize(rfftp_plan plan)
        {
            int twsize = 0, l1 = 1;
            for (int k = 0; k < plan.nfct; ++k)
            {
                int ip = plan.fct[k].fct, ido = plan.length / (l1 * ip);
                twsize += (ip - 1) * (ido - 1);
                if (ip > 5) twsize += 2 * ip;
                l1 *= ip;
            }

            return twsize;
        }

        static void rfftp_comp_twiddle(rfftp_plan plan)
        {
            int length = plan.length;
            double[] twid = new double[2 * length];
            sincos_2pibyn_half(length, twid);
            int l1 = 1;
            int ptr = 0;
            for (int k = 0; k < plan.nfct; ++k)
            {
                int ip = plan.fct[k].fct, ido = length / (l1 * ip);
                if (k < plan.nfct - 1) // last factor doesn't need twiddles
                {
                    plan.fct[k].tw = ptr;
                    Span<double> tw = plan.mem.AsSpan(ptr);
                    ptr += (ip - 1) * (ido - 1);
                    for (int j = 1; j < ip; ++j)
                    {
                        for (int i = 1; i <= (ido - 1) / 2; ++i)
                        {
                            tw[(j - 1) * (ido - 1) + 2 * i - 2] = twid[2 * j * l1 * i];
                            tw[(j - 1) * (ido - 1) + 2 * i - 1] = twid[2 * j * l1 * i + 1];
                        }
                    }
                }

                if (ip > 5) // special factors required by *g functions
                {
                    plan.fct[k].tws = ptr;
                    Span<double> tws = plan.mem.AsSpan(ptr);
                    ptr += 2 * ip;
                    tws[0] = 1.0;
                    tws[1] = 0.0;
                    for (int i = 1; i <= (ip >> 1); ++i)
                    {
                        tws[2 * i] = twid[2 * i * (length / ip)];
                        tws[2 * i + 1] = twid[2 * i * (length / ip) + 1];
                        tws[2 * (ip - i)] = twid[2 * i * (length / ip)];
                        tws[2 * (ip - i) + 1] = -twid[2 * i * (length / ip) + 1];
                    }
                }

                l1 *= ip;
            }

            //DEALLOC(twid);
        }

        static rfftp_plan make_rfftp_plan(int length)
        {
            if (length == 0) return null;
            rfftp_plan plan = new rfftp_plan();

            plan.length = length;
            plan.nfct = 0;
            plan.mem = null;
#if !NET9_0_OR_GREATER
            if (plan.fct == null) // if it's not an inline struct array...
            {
                plan.fct = new rfftp_fctdata[NFCT];
            }
#endif

            for (int i = 0; i < NFCT; ++i)
            {
                plan.fct[i] = new rfftp_fctdata(0, 0, 0);
            }

            if (length == 1)
            {
                return plan;
            }

            if (!rfftp_factorize(plan))
            {
                return null;
            }

            int tws = rfftp_twsize(plan);
            plan.mem = new double[tws];
            rfftp_comp_twiddle(plan);

            //if (NAIL_TEST) Console.WriteLine("rfftp plan mem");
            //if (NAIL_TEST) NailTestPrintDoubleArray(plan.mem.AsSpan(0, tws));

            return plan;
        }

        static void destroy_rfftp_plan(rfftp_plan plan)
        {
        }

        public class fftblue_plan
        {
            public int n;
            public int n2;
            public cfftp_plan plan;
            public double[] mem;
            public int bk; // indexes into mem
            public int bkf; // indexes into mem

            public IEnumerable<int> Factors()
            {
                return plan.Factors();
            }
        }

        static fftblue_plan make_fftblue_plan(int length)
        {
            fftblue_plan plan = new fftblue_plan();
            plan.n = length;
            plan.n2 = good_size(plan.n * 2 - 1);
            plan.mem = new double[2 * plan.n + 2 * plan.n2];
            plan.bk = 0;
            plan.bkf = plan.bk + 2 * plan.n;

            /* initialize b_k */
            double[] tmp = new double[4 * plan.n];
            sincos_2pibyn(2 * plan.n, tmp);
            Span<double> bk = plan.mem.AsSpan(plan.bk);
            Span<double> bkf = plan.mem.AsSpan(plan.bkf);
            bk[0] = 1;
            bk[1] = 0;

            int coeff = 0;
            for (int m = 1; m < plan.n; ++m)
            {
                coeff += 2 * m - 1;
                if (coeff >= 2 * plan.n) coeff -= 2 * plan.n;
                bk[2 * m] = tmp[2 * coeff];
                bk[2 * m + 1] = tmp[2 * coeff + 1];
            }

            /* initialize the zero-padded, Fourier transformed b_k. Add normalisation. */
            double xn2 = 1.0 / plan.n2;
            bkf[0] = bk[0] * xn2;
            bkf[1] = bk[1] * xn2;
            for (int m = 2; m < 2 * plan.n; m += 2)
            {
                bkf[m] = bkf[2 * plan.n2 - m] = bk[m] * xn2;
                bkf[m + 1] = bkf[2 * plan.n2 - m + 1] = bk[m + 1] * xn2;
            }

            for (int m = 2 * plan.n; m <= (2 * plan.n2 - 2 * plan.n + 1); ++m)
            {
                bkf[m] = 0.0;
            }

            plan.plan = make_cfftp_plan(plan.n2);
            if (plan.plan == null)
            {
                return null;
            }

            cfftp_forward(plan.plan, bkf, 1.0);

            //DEALLOC(tmp);
            return plan;
        }

        static void destroy_fftblue_plan(fftblue_plan plan)
        {
        }

        static void fftblue_fft(fftblue_plan plan, Span<double> c, int isign, double fct)
        {
            int n = plan.n;
            int n2 = plan.n2;
            Span<double> bk = plan.mem.AsSpan(plan.bk);
            Span<double> bkf = plan.mem.AsSpan(plan.bkf);
            Span<double> akf = new double[2 * n2];

            /* initialize a_k and FFT it */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
                    akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
                }
            }
            else
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
                    akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
                }
            }
            for (int m = 2 * n; m < 2 * n2; ++m)
            {
                akf[m] = 0;
            }

            cfftp_forward(plan.plan, akf, fct);

            /* do the convolution */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n2; m += 2)
                {
                    double im = -akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
                    akf[m] = akf[m] * bkf[m] + akf[m + 1] * bkf[m + 1];
                    akf[m + 1] = im;
                }
            }
            else
            {
                for (int m = 0; m < 2 * n2; m += 2)
                {
                    double im = akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
                    akf[m] = akf[m] * bkf[m] - akf[m + 1] * bkf[m + 1];
                    akf[m + 1] = im;
                }
            }

            /* inverse FFT */
            cfftp_backward(plan.plan, akf, 1.0);

            /* multiply by b_k */
            if (isign > 0)
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
                    c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
                }
            }
            else
            {
                for (int m = 0; m < 2 * n; m += 2)
                {
                    c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
                    c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
                }
            }

            //DEALLOC(akf);
        }

        static void cfftblue_backward(fftblue_plan plan, Span<double> c, double fct)
        {
            fftblue_fft(plan, c, 1, fct);
        }

        static void cfftblue_forward(fftblue_plan plan, Span<double> c, double fct)
        {
            fftblue_fft(plan, c, -1, fct);
        }

        static void rfftblue_backward(fftblue_plan plan, Span<double> c, double fct)
        {
            int n = plan.n;
            double[] tmp = new double[2 * n];
            tmp[0] = c[0];
            tmp[1] = 0.0;
            c.Slice(1, (n - 1)).CopyTo(tmp.AsSpan(2));
            //memcpy(tmp + 2, c + 1, (n - 1) * sizeof(double));

            if ((n & 1) == 0) tmp[n + 1] = 0.0;
            for (int m = 2; m < n; m += 2)
            {
                tmp[2 * n - m] = tmp[m];
                tmp[2 * n - m + 1] = -tmp[m + 1];
            }

            fftblue_fft(plan, tmp, 1, fct);
            for (int m = 0; m < n; ++m)
            {
                c[m] = tmp[2 * m];
            }

            //DEALLOC(tmp);
        }

        static void rfftblue_forward(fftblue_plan plan, Span<double> c, double fct)
        {
            int n = plan.n;
            double[] tmp = new double[2 * n];
            for (int m = 0; m < n; ++m)
            {
                tmp[2 * m] = c[m];
                tmp[2 * m + 1] = 0.0;
            }

            fftblue_fft(plan, tmp, -1, fct);

            c[0] = tmp[0];
            tmp.AsSpan(2, (n - 1)).CopyTo(c.Slice(1));
            //memcpy(c + 1, tmp + 2, (n - 1) * sizeof(double));

            //DEALLOC(tmp);
        }

        public class cfft_plan
        {
            public cfftp_plan packplan;
            public fftblue_plan blueplan;

            public IEnumerable<int> Factors()
            {
                if (packplan != null)
                {
                    return packplan.Factors();
                }
                else
                {
                    return blueplan.Factors();
                }
            }
        }

        public static cfft_plan make_cfft_plan(int length)
        {
            if (length == 0) return null;
            cfft_plan plan = new cfft_plan();

            //plan.blueplan = null;
            //plan.packplan = null;

            if ((length < 50) || (largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                plan.packplan = make_cfftp_plan(length);
                if (plan.packplan == null)
                {
                    return null;
                }

                return plan;
            }

            double comp1 = cost_guess(length);
            double comp2 = 2 * cost_guess(good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */

            if (comp2 < comp1) // use Bluestein
            {
                plan.blueplan = make_fftblue_plan(length);

                if (plan.blueplan == null)
                {
                    return null;
                }
            }
            else
            {
                plan.packplan = make_cfftp_plan(length);

                if (plan.packplan == null)
                {
                    return null;
                }
            }

            return plan;
        }

        public static void destroy_cfft_plan(cfft_plan plan)
        {
        }

        public static void cfft_backward(cfft_plan plan, Span<double> c, double fct)
        {
            if (plan.packplan != null)
            {
                cfftp_backward(plan.packplan, c, fct);
            }
            else
            {
                cfftblue_backward(plan.blueplan, c, fct);
            }
        }

        public static void cfft_forward(cfft_plan plan, Span<double> c, double fct)
        {
            if (plan.packplan != null)
            {
                cfftp_forward(plan.packplan, c, fct);
            }
            else
            {
                cfftblue_forward(plan.blueplan, c, fct);
            }
        }

        public class rfft_plan
        {
            public rfftp_plan packplan;
            public fftblue_plan blueplan;
            public IEnumerable<int> Factors()
            {
                if (packplan != null)
                {
                    return packplan.Factors();
                }
                else
                {
                    return blueplan.Factors();
                }
            }
        }

        public static rfft_plan make_rfft_plan(int length)
        {
            if (length == 0) return null;
            rfft_plan plan = new rfft_plan();
            //plan.blueplan = null;
            //plan.packplan = null;

            if ((length < 50) || (largest_prime_factor(length) <= Math.Sqrt(length)))
            {
                plan.packplan = make_rfftp_plan(length);
                if (plan.packplan == null)
                {
                    return null;
                }

                return plan;
            }

            double comp1 = 0.5 * cost_guess(length);
            double comp2 = 2 * cost_guess(good_size(2 * length - 1));
            comp2 *= 1.5; /* fudge factor that appears to give good overall performance */
            if (comp2 < comp1) // use Bluestein
            {
                plan.blueplan = make_fftblue_plan(length);
                if (plan.blueplan == null)
                {
                    return null;
                }
            }
            else
            {
                plan.packplan = make_rfftp_plan(length);
                if (plan.packplan == null)
                {
                    return null;
                }
            }

            return plan;
        }

        public static void destroy_rfft_plan(rfft_plan plan)
        {
        }

        public static int rfft_length(rfft_plan plan)
        {
            if (plan.packplan != null)
            {
                return plan.packplan.length;
            }

            return plan.blueplan.n;
        }

        public static int cfft_length(cfft_plan plan)
        {
            if (plan.packplan != null)
            {
                return plan.packplan.length;
            }

            return plan.blueplan.n;
        }

        public static void rfft_backward(rfft_plan plan, Span<double> c, double fct)
        {
            if (plan.packplan != null)
            {
                rfftp_backward(plan.packplan, c, fct);
            }

            else // if (plan.blueplan)
                rfftblue_backward(plan.blueplan, c, fct);
        }

        public static void rfft_forward(rfft_plan plan, Span<double> c, double fct)
        {
            if (plan.packplan != null)
            {
                rfftp_forward(plan.packplan, c, fct);
            }
            else // if (plan.blueplan)
            {
                rfftblue_forward(plan.blueplan, c, fct);
            }
        }
    }
}