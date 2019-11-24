using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Weinstein_and_Fiker_methods
{
    class WeinsteinMethod
    {
        private double _left;
        private double _right;
        private int _n;
        private double _step;

        private Func<double, double> _function;

        private double[] _disceteFunction;
        private double[,] _T;
        private double[,] _A;
        private double[,] _B;
        private double[,] _C;

        public WeinsteinMethod(Func<double, double> function, int n, double left, double right)
        {
            _function = function;
            _n = n;
            _left = left;
            _right = right;
            _step = (_right - _left) / (_n - 1);

            // Calculating discrete functon
            _disceteFunction = new double[_n];
            for (int i = 0; i < _n; ++i)
            {
                _disceteFunction[i] = _function(_left + i * _step);
            }

            // Creating matrix T
            _T = new double[_n, _n];
            for (int i = 0; i < _n; ++i)
            {
                _T[i, i] = -2;
                if (i + 1 < _n)
                {
                    _T[i, i + 1] = 1;
                    _T[i + 1, i] = 1;
                }
            }

            _A = CreateMatrixA();
            _B = CreateMatrixB();
            _C = CreateMatrixC();
        }

        private double[,] CreateMatrixA()
        {
            double[,] result = new double[_n, _n];
            for (int i = 0; i < _n; ++i)
            {
                result[i, i] = _disceteFunction[i];
            }
            return _T.Mult(result);
        }

        private double[,] CreateMatrixB()
        {
            double[,] result = CreateMatrixA();
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    result[i, j] -= 0.05;
                }
            }
            return result;
        }

        private double[,] CreateMatrixC()
        {
            double[,] result = new double [_n, _n];
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    result[i, j] = _A[i, j] - _B[i, j];
                }
            }
            return result;
        }

        private double[] CreateFkArray(int k)
        {
            var fk = new double[_n];
            if (k == 0)
            {
                for (int i = 0; i < _n; ++i)
                {
                    fk[i] = 1;
                }
            }
            else
            {
                for (int i = 0; i < _n; ++i)
                {
                    fk[i] = Math.Sin(_left + i * _step);
                }
            }
            return fk;
        }

        private double[] CreateGkArray(double[] fk)
        {
            // Matrix * vector
            return Matrix.Dot(_C, fk);
        }

        private double[] GetCalculatedEigenValues(double[,] A, int m, double[,] Am)
        {
            var inverse = new Accord.Math.Decompositions.EigenvalueDecomposition(A.Inverse(), false, true);
            var eigVals = inverse.RealEigenvalues;
            double dev = ((double)m).GetDev();
            for (int i = 0; i < eigVals.Length; ++i)
            {
                eigVals[i] = 1 / eigVals[i];
                eigVals[i] = i % 2 != 0 ? eigVals[i] + dev : eigVals[i] - dev;
            }
            return eigVals.OrderBy(c => c).ToList().ToArray();
        }

        private double[,] GetCmOperator(double[,] C, int m)
        {
            List<double[]> fkArray = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                fkArray.Add(CreateFkArray(i));
            }

            var gkArray = new List<double[]>();
            for (int i = 0; i < m; ++i)
            {
                gkArray.Add(CreateGkArray(fkArray[i]));
            }

            var bij = new double[m, m];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    double sum = 0;
                    for (int k = 0; k < fkArray[j].Length; ++k)
                    {
                        sum += fkArray[j][k] * gkArray[i][k];
                    }
                    bij[i, j] += sum;
                }
            }

            var Cm = new double[_n, _n];
            for (int i = 0; i < m; ++i)
            {
                var gki = new double[_n];
                for (int k = 0; k < _n; ++k)
                {
                    gki[k] = gkArray[i][k];
                }

                for (int j = 0; j < m; ++j)
                {
                    for (int z = 0; z < _n; ++z)
                    {
                        gki[z] = gki[z] * bij[i, j];
                    }

                    double[,] matr = new double[_n, _n];
                    for (int l = 0; l < _n; ++l)
                    {
                        for (int o = 0; o < _n; ++o)
                        {
                            matr[l, o] = gki[l] * gkArray[j][o];
                        }
                    }

                    Cm = matr;
                }
            }
            return Cm;
        }

        public void Calculate(int m, out double[] realValues, out double[] calculatedValues)
        {
            var Cm = GetCmOperator(_C, m);
            var Am = new double[_n, _n];
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    Am[i, j] = _B[i, j] + Cm[i, j];
                }
            }
            calculatedValues = GetCalculatedEigenValues(_A, m, Am).OrderBy(c => c).ToArray();

            var eigReal = new Accord.Math.Decompositions.EigenvalueDecomposition(_A, false, true);
            realValues = eigReal.RealEigenvalues.OrderBy(c => c).ToArray();
        }
    }
}
