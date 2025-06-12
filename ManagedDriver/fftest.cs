using PocketFFT;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ManagedDriver
{
    public static class fftest
    {
        private const int maxlen = 8196;

        private static void fill_random(Span<double> data, int length)
        {
            Random rand = new Random(77143423);
            for (int m = 0; m < length; ++m)
            {
                data[m] = rand.NextDouble() - 0.5;
            }
            //data[0] = 1.0;
        }

        private static double errcalc(Span<double> data, Span<double> odata, int length)
        {
            double sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m] - odata[m]) * (data[m] - odata[m]);
                sum += odata[m] * odata[m];
            }

            return Math.Sqrt(errsum / sum);
        }

        public static int test_real()
        {
            double[] data = new double[maxlen];
            double[] odata = new double[maxlen];
            const double epsilon = 2e-15;
            int ret = 0;
            fill_random(odata, maxlen);
            double errsum = 0;
            for (int length = 1; length <= maxlen; ++length)
            {
                //Console.WriteLine("testing real " + length);
                odata.AsSpan(0, length).CopyTo(data);
                PocketFFTBaseline.rfft_plan plan = PocketFFTBaseline.make_rfft_plan(length);
                PocketFFTBaseline.rfft_forward(plan, data, 1.0);
                PocketFFTBaseline.rfft_backward(plan, data, 1.0 / length);
                PocketFFTBaseline.destroy_rfft_plan(plan);
                double err = errcalc(data, odata, length);
                if (err > epsilon)
                {
                    Console.WriteLine($"problem at real length {length}: {err} - factors {string.Join(',', plan.Factors())}");
                    ret = 1;
                }
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
            return ret;
        }

        public static int test_complex()
        {
            double[] data = new double[2 * maxlen];
            double[] odata = new double[2 * maxlen];
            fill_random(odata, 2 * maxlen);
            const double epsilon = 2e-15;
            int ret = 0;
            double errsum = 0;

            for (int length = 1; length <= maxlen; ++length)
            {
                //Console.WriteLine("testing cmplx " + length);
                odata.AsSpan(0, 2 * length).CopyTo(data);
                PocketFFTBaseline.cfft_plan plan = PocketFFTBaseline.make_cfft_plan(length);
                PocketFFTBaseline.cfft_forward(plan, data, 1.0);
                PocketFFTBaseline.cfft_backward(plan, data, 1.0 / length);
                PocketFFTBaseline.destroy_cfft_plan(plan);
                double err = errcalc(data, odata, 2 * length);
                if (err > epsilon)
                {
                    Console.WriteLine($"problem at complex length {length}: {err} - factors {string.Join(',', plan.Factors())}");
                    ret = 1;
                }
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
            return ret;
        }

        public static void TestReal(int length)
        {
            double[] data = new double[length];
            double[] odata = new double[length];
            fill_random(odata, length);
            odata.AsSpan(0, length).CopyTo(data);
            PocketFFTBaseline.rfft_plan plan = PocketFFTBaseline.make_rfft_plan(length);
            PocketFFTBaseline.rfft_forward(plan, data, 1.0);
            PocketFFTBaseline.rfft_backward(plan, data, 1.0 / length);
            PocketFFTBaseline.destroy_rfft_plan(plan);
            double err = errcalc(data, odata, length);
            Console.WriteLine($"err: {err}");
        }
    }
}
