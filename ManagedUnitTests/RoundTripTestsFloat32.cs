using PocketFFT;

namespace ManagedUnitTests
{
    [TestClass]
    public class RoundTripTestsFloat32   
    {
        private const int MAX_LENGTH = 8192;
        private const float epsilon = 2e-6f;

        [TestMethod]
        public void TestReal()
        {
            float errsum = 0;
            for (int length = 1; length <= MAX_LENGTH; ++length)
            {
                //Console.WriteLine("testing real " + length);
                float[] data = new float[2 * length];
                float[] odata = new float[2 * length];
                fill_random(odata, length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IReal1DFFTPlanFloat32 plan = FFTPlanFactory.Create1DRealFFTPlanFloat32(length))
                {
                    plan.Forward(data, 1.0f);
                    plan.Backward(data, 1.0f / length);
                }

                float err = errcalc(data, odata, length);
                Assert.IsTrue(err < epsilon, $"problem at real length {length}: {err}");
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
        }

        [TestMethod]
        public void TestComplex()
        {
            float errsum = 0;

            for (int length = 1; length <= MAX_LENGTH; ++length)
            {
                //Console.WriteLine("testing cmplxF " + length);
                cmplxF[] data = new cmplxF[length];
                cmplxF[] odata = new cmplxF[length];
                fill_random(odata, length);
                odata.AsSpan(0, length).CopyTo(data);
                using (IComplex1DFFTPlanFloat32 plan = FFTPlanFactory.Create1DComplexFFTPlanFloat32(length))
                {
                    plan.Forward(data, 1.0f);
                    plan.Backward(data, 1.0f / length);
                }
                float err = errcalc(data, odata, length);
                Assert.IsTrue(err < epsilon, $"problem at complex length {length}: {err}");
                errsum += err;
            }

            Console.WriteLine($"errsum: {errsum}");
        }

        private static void fill_random(Span<float> data, int length)
        {
            for (int m = 0; m < length; ++m)
            {
                data[m] = (float)(Random.Shared.NextDouble() - 0.5);
            }
        }

        private static void fill_random(Span<cmplxF> data, int length)
        {
            for (int m = 0; m < length; ++m)
            {
                data[m].r = (float)(Random.Shared.NextDouble() - 0.5);
                data[m].i = (float)(Random.Shared.NextDouble() - 0.5);
            }
        }

        private static float errcalc(Span<float> data, Span<float> odata, int length)
        {
            float sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m] - odata[m]) * (data[m] - odata[m]);
                sum += odata[m] * odata[m];
            }

            return (float)Math.Sqrt(errsum / sum);
        }

        private static float errcalc(Span<cmplxF> data, Span<cmplxF> odata, int length)
        {
            float sum = 0, errsum = 0;
            for (int m = 0; m < length; ++m)
            {
                errsum += (data[m].r - odata[m].r) * (data[m].r - odata[m].r);
                errsum += (data[m].i - odata[m].i) * (data[m].i - odata[m].i);
                sum += odata[m].r * odata[m].r;
                sum += odata[m].i * odata[m].i;
            }

            return (float)Math.Sqrt(errsum / sum);
        }
    }
}