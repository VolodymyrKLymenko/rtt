using System;

namespace ExactDerivativesAndInitApproach
{
    class Program
    {
        public static double F(double x)
        {
            return -((4.0 - x * x) * (4.0 - x * x));
        }

        static void Main(string[] args)
        {
            ExactDerivatives exactD = new ExactDerivatives(F, -1, 1, 50);
            exactD.Calculate(0.07, 0.0001, out double calculated, out int numOfIterations, out double[] exactValues);

            Console.WriteLine("Calculating exact derivatives:");
            Console.WriteLine("Exact values:         " + string.Join("\t", exactValues));
            Console.WriteLine("Calculated value:     " + calculated);
            Console.WriteLine("Number of iterations: " + numOfIterations);

            Console.WriteLine("__________________________________________________________\n");

            /*InitialApproximation appx = new InitialApproximation(F, 0, 1, 50);
            appx.Calculate(0, 0.2, out int numberOfLambdas, out double[] lambdas);
            Console.WriteLine("Looking for initial approximation: found " + numberOfLambdas + ":");

            for (int i = 0; i < numberOfLambdas; ++i)
            {
                Console.WriteLine(lambdas[i].ToString());
            }*/
            
            Console.ReadKey();
        }
    }
}
