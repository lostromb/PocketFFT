using PocketFFT;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ManagedUnitTests
{
    [TestClass]
    public class IntrinsicTests
    {
        [TestMethod]
        public void TestSpanRefEquals()
        {
            byte[] array = new byte[1000];
            byte[] array2 = new byte[1000];
            Span<byte> a = array.AsSpan();
            Span<byte> b = array.AsSpan();
            Assert.IsTrue(Intrinsics.SpanRefEquals(a, b));
            Assert.IsFalse(Intrinsics.SpanRefEquals(a, b.Slice(1)));
            Assert.IsTrue(Intrinsics.SpanRefEquals(a.Slice(1), b.Slice(1)));
            Assert.IsTrue(Intrinsics.SpanRefEquals(array.AsSpan(), array.AsSpan()));
            Assert.IsFalse(Intrinsics.SpanRefEquals(array.AsSpan(), array2.AsSpan()));
            Assert.IsFalse(Intrinsics.SpanRefEquals(a, Span<byte>.Empty));
            Assert.IsFalse(Intrinsics.SpanRefEquals(Span<byte>.Empty, b));
            Assert.IsTrue(Intrinsics.SpanRefEquals(Span<byte>.Empty, Span<byte>.Empty));
        }

        [TestMethod]
        public void TestCastSingleToDouble()
        {
            for (int length = 0; length < 1000; length++)
            {
                float[] input = new float[length];
                double[] output = new double[length];
                for (int c = 0; c < length; c++)
                {
                    input[c] = Random.Shared.NextSingle();
                }

                Intrinsics.CastSingleToDouble(input.AsSpan(), output.AsSpan());

                for (int c = 0; c < length; c++)
                {
                    Assert.AreEqual(input[c], (float)output[c], 0.00001f);
                }
            }
        }

        [TestMethod]
        public void TestCastDoubleToSingle()
        {
            for (int length = 0; length < 1000; length++)
            {
                double[] input = new double[length];
                float[] output = new float[length];
                for (int c = 0; c < length; c++)
                {
                    input[c] = Random.Shared.NextDouble();
                }

                Intrinsics.CastDoubleToSingle(input.AsSpan(), output.AsSpan());

                for (int c = 0; c < length; c++)
                {
                    Assert.AreEqual(input[c], (double)output[c], 0.00001f);
                }
            }
        }
    }
}
