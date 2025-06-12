using PocketFFT;

namespace ManagedUnitTests
{
    [TestClass]
    public class RoundTripTestsFloat64
    {
        private const int MAX_LENGTH = 8192;
        private const double epsilon = 2e-15;

        [TestMethod]
        public void TestReal()
        {
            double errsum = 0;
            for (int length = 1; length <= MAX_LENGTH; ++length)
            {
                //Console.WriteLine("testing real " + length);
                double[] data = new double[2 * length];
                double[] odata = new double[2 * length];
                fill_random(odata, length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IReal1DFFTPlanFloat64 plan = FFTPlanFactory.Create1DRealFFTPlanFloat64(length))
                {
                    plan.Forward(data, 1.0);
                    plan.Backward(data, 1.0 / length);
                }

                double err = errcalc(data, odata, length);
                Assert.IsTrue(err < epsilon, $"problem at real length {length}: {err}");
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
        }

        [TestMethod]
        public void TestComplex()
        {
            double errsum = 0;

            for (int length = 1; length <= MAX_LENGTH; ++length)
            {
                //Console.WriteLine("testing cmplx " + length);
                cmplx[] data = new cmplx[length];
                cmplx[] odata = new cmplx[length];
                fill_random(odata, length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IComplex1DFFTPlanFloat64 plan = FFTPlanFactory.Create1DComplexFFTPlanFloat64(length))
                {
                    plan.Forward(data, 1.0);
                    plan.Backward(data, 1.0 / length);
                }
                double err = errcalc(data, odata, length);
                Assert.IsTrue(err < epsilon, $"problem at complex length {length}: {err}");
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
        }

        private static void fill_random(Span<double> data, int length)
        {
            for (int m = 0; m < length; ++m)
            {
                data[m] = Random.Shared.NextDouble() - 0.5;
            }
        }

        private static void fill_random(Span<cmplx> data, int length)
        {
            for (int m = 0; m < length; ++m)
            {
                data[m].r = Random.Shared.NextDouble() - 0.5;
                data[m].i = Random.Shared.NextDouble() - 0.5;
            }
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

        private static double errcalc(Span<cmplx> data, Span<cmplx> odata, int length)
        {
            double sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m].r - odata[m].r) * (data[m].r - odata[m].r);
                errsum += (data[m].i - odata[m].i) * (data[m].i - odata[m].i);
                sum += odata[m].r * odata[m].r;
                sum += odata[m].i * odata[m].i;
            }

            return Math.Sqrt(errsum / sum);
        }
    }
}