using BenchmarkDotNet.Running;
using ManagedDriver;
using PocketFFT;

namespace Driver
{
    internal static class Program
    {
        public static void Main(string[] args)
        {
            //fftest.TestReal(189);
            //fftest.TestReal(190);
            //fftest.TestReal(191);
            //fftest.test_real();
            //fftest.test_complex();
            BenchmarkRunner.Run<Benchmarks>();
        }
    }
}
