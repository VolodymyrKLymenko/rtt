using System;
using System.Linq;

namespace ExactDerivativesAndInitApproach
{
    class ExactDerivatives
    {
        private int _n;
        private double _h;

        private double[] _discreteFunction;
        private double[,] _T;
        private double[,] _A;

        public ExactDerivatives(Func<double, double> function, double left, double right, int n)
        {
            _n = n;
            _h = (right - left) / (_n - 1);

            // Crating the dicrete function
            _discreteFunction = new double[_n];
            for (int i = 0; i < _n; ++i)
            {
                _discreteFunction[i] = function(left + i * _h);
            }

            // Creating T
            _T = new double[n, n];
            for (int i = 0; i < _n; ++i)
            {
                _T[i, i] = -2;
                if (i + 1 < _n)
                {
                    _T[i, i + 1] = 1;
                    _T[i + 1, i] = 1;
                }
            }

            // Creating A
            double[,] discreteFMatrix = new double[_n, _n];
            for (int i = 0; i < _n; ++i)
            {
                discreteFMatrix[i, i] = _discreteFunction[i];
            }
            _A = _T.Multiply(discreteFMatrix);
        }

        private double[,] CreateMatrixD(double lambda)
        {
            double[,] diagonalMatrix = new double[_n, _n];
            diagonalMatrix = diagonalMatrix.SetDiagonalElements(lambda);
            return _A.Subtract(diagonalMatrix);
        }

        private static double GetDeltaLambda(double[,] d, int n)
        {
            double[,] b = new double[n, n];
            b.SetDiagonalElements(1);

            double[,] u = new double[n, n];
            double[,] l = new double[n, n];
            double[,] v = new double[n, n];
            double[,] m = new double[n, n];
            
            for (int r = 0; r < n; ++r)
            {
                int k = r;
                for (int i = r + 1; i < n; ++i, ++k)
                {
                    u[r, k] =  d[r, k] - l.GetRow(r).Multiply(u.GetColumn(k)).SumFirstElements(r);

                    l[i, r] = (d[i, r] - l.GetRow(i).Multiply(u.GetColumn(r)).SumFirstElements(r)) / u[r, r];

                    v[r, k] =  b[r, k] - m.GetRow(r).Multiply(u.GetColumn(k))
                        .Sum(l.GetRow(r).Multiply(v.GetColumn(k))).SumFirstElements(r);

                    m[i, r] = (b[i, r] - m.GetRow(i).Multiply(u.GetColumn(r))
                        .Sum(l.GetRow(i).Multiply(v.GetColumn(r))).SumFirstElements(r) - l[i, r] * v[r, r]) / u[r, r];
                }
                k = n - 1;

                u[r, k] = d[r, k] - l.GetRow(r).Multiply(u.GetColumn(k))
                    .SumFirstElements(r);
                v[r, k] = b[r, k] - m.GetRow(r).Multiply(u.GetColumn(k))
                    .Sum(l.GetRow(r).Multiply(v.GetColumn(k))).SumFirstElements(r);
            }

            double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                sum += v[i, i] / u[i, i];
            }

            return 1 / sum;
        }

        public void Calculate(
            double initialApproach,
            double deviation,
            out double calculatedValue,
            out int numberOfIterations,
            out double[] exactValues)
        {
           
            double lambdaPrev = initialApproach;
            double lambdaNext = 1.3;
            double deltaLambda = double.MaxValue;
            numberOfIterations = 0;

            while (Math.Abs(deltaLambda) > deviation)
            {
                double[,] D = CreateMatrixD(lambdaPrev);
                deltaLambda = GetDeltaLambda(D, _n);

                lambdaNext = lambdaPrev + deltaLambda;
                lambdaPrev = lambdaNext;
                ++numberOfIterations;
            }

            calculatedValue = lambdaNext;

            var decA = new Accord.Math.Decompositions.EigenvalueDecomposition(_A, false, true);
            exactValues = decA.RealEigenvalues.OrderBy(c => c).ToArray();
        }
    }
}
