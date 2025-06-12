using BenchmarkDotNet.Attributes;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using static Driver.Benchmarks;

namespace Driver
{
    [DisassemblyDiagnoser]
    public class Benchmarks
    {
        private Complex[]? cc;
        private Complex[]? ch;
        private Complex[]? wa;

        [GlobalSetup]
        public void Setup()
        {
            cc = new Complex[65536];
            ch = new Complex[65536];
            wa = new Complex[65536];
            Span<Complex> span = cc;
            foreach (ref Complex vec in span)
            {
                vec.r = (float)Random.Shared.NextDouble();
                vec.i = (float)Random.Shared.NextDouble();
            }

            span = ch;
            foreach (ref Complex vec in span)
            {
                vec.r = (float)Random.Shared.NextDouble();
                vec.i = (float)Random.Shared.NextDouble();
            }

            span = wa;
            foreach (ref Complex vec in span)
            {
                vec.r = (float)Random.Shared.NextDouble();
                vec.i = (float)Random.Shared.NextDouble();
            }
        }

        [Benchmark]
        public void Pass2Inline()
        {
            pass2b(2, 1024, cc.AsSpan(), ch.AsSpan(), wa.AsSpan());
        }

        [Benchmark]
        public void Pass2RefFunc()
        {
            pass2b_reffunc(2, 1024, cc.AsSpan(), ch.AsSpan(), wa.AsSpan());
        }

        [Benchmark]
        public void Pass2Assign()
        {
            pass2b_assign(2, 1024, cc.AsSpan(), ch.AsSpan(), wa.AsSpan());
        }

        internal struct Complex
        {
            internal double r;
            internal double i;

            public Complex(double r, double i)
            {
                this.r = r;
                this.i = i;
            }
        }

        private static void PMC(ref Complex a, ref Complex b, ref Complex c, ref Complex d)
        {
            a.r = c.r + d.r;
            a.i = c.i + d.i;
            b.r = c.r - d.r;
            b.i = c.i - d.i;
        }

        private static void A_EQ_B_MUL_C(ref Complex a, ref Complex b, ref Complex c)
        {
            a.r = b.r * c.r - b.i * c.i;
            a.i = b.r * c.i + b.i * c.r;
        }
        private static void pass2b_reffunc(int ido, int l1, Span<Complex> cc, Span<Complex> ch, Span<Complex> wa)
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
                        Complex t = default;
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

        private static void pass2b_assign(int ido, int l1, Span<Complex> cc, Span<Complex> ch, Span<Complex> wa)
        {
            const int cdim = 2;

            if (ido == 1)
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))] = new Complex(
                        cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r,
                        cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i);
                    ch[(0) + ido * ((k) + l1 * (1))] = new Complex(
                        cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r,
                        ch[(0) + ido * ((k) + l1 * (1))].i = cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i);
                }
            }
            else
            {
                for (int k = 0; k < l1; ++k)
                {
                    ch[(0) + ido * ((k) + l1 * (0))] = new Complex(
                        cc[(0) + ido * ((0) + cdim * (k))].r + cc[(0) + ido * ((1) + cdim * (k))].r,
                        cc[(0) + ido * ((0) + cdim * (k))].i + cc[(0) + ido * ((1) + cdim * (k))].i);
                    ch[(0) + ido * ((k) + l1 * (1))] = new Complex(
                        cc[(0) + ido * ((0) + cdim * (k))].r - cc[(0) + ido * ((1) + cdim * (k))].r,
                        cc[(0) + ido * ((0) + cdim * (k))].i - cc[(0) + ido * ((1) + cdim * (k))].i);

                    for (int i = 1; i < ido; ++i)
                    {
                        Complex t;
                        ch[(i) + ido * ((k) + l1 * (0))] = new Complex(
                            cc[(i) + ido * ((0) + cdim * (k))].r + cc[(i) + ido * ((1) + cdim * (k))].r,
                            cc[(i) + ido * ((0) + cdim * (k))].i + cc[(i) + ido * ((1) + cdim * (k))].i);
                        t = new Complex(
                            cc[(i) + ido * ((0) + cdim * (k))].r - cc[(i) + ido * ((1) + cdim * (k))].r,
                            cc[(i) + ido * ((0) + cdim * (k))].i - cc[(i) + ido * ((1) + cdim * (k))].i);
                        ch[(i) + ido * ((k) + l1 * (1))] = new Complex(
                            wa[(i) - 1 + (0) * (ido - 1)].r * t.r - wa[(i) - 1 + (0) * (ido - 1)].i * t.i,
                            wa[(i) - 1 + (0) * (ido - 1)].r * t.i + wa[(i) - 1 + (0) * (ido - 1)].i * t.r);
                    }
                }
            }
        }

        private static void pass2b(int ido, int l1, Span<Complex> cc, Span<Complex> ch, Span<Complex> wa)
        {
            const int cdim = 2;

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
                        Complex t;
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
    }
}
