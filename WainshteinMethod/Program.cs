using System;
using System.Linq;

namespace WainshteinMethod
{
    public class Program
    {
        static void Main(string[] args)
        {
            var left = -1;
            var right = 1;
            var discretizationCount = 100;
            var F = new Func<double, double>(x => -((4.0 - x * x) * (4.0 - x * x)));

            var wainshteingExecution = new WainshteinMethod(F, discretizationCount, left, right);

            wainshteingExecution.Calculate(10, out var realValues, out var calculatedValues);

            Console.WriteLine("Real values");
            foreach (var realValue in realValues.OrderByDescending(v => v))
            {
                Console.Write($"{realValue}, ");
            }

            Console.WriteLine("\n\n\nCalculated values");
            foreach (var calculatedValue in calculatedValues.OrderByDescending(v => v))
            {
                Console.Write($"{calculatedValue}, ");
            }

            Console.ReadLine();
        }
    }
}
