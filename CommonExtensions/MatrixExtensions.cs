using MathNet.Numerics.LinearAlgebra;

namespace ExactDerivativesAndInitApproach
{
    public static class Extensions
    {
        public static Vector<double> Multiply(this Vector<double> first, Vector<double> second, double constant = 1.0)
        {
            var result = Vector<double>.Build.Dense(first.Count);

            for (int i = 0; i < first.Count; ++i)
            {
                result[i] = first[i] * second[i] * constant;
            }
            return result;
        }

        public static Matrix<double> Multiply(this Matrix<double> first, Matrix<double> second)
        {
            int rows = first.RowCount;
            int cols = first.ColumnCount;
            var result = Matrix<double>.Build.Dense(rows, cols);

            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    for (int k = 0; k < cols; ++k)
                    {
                        result[i, j] += first[i, k] * second[k, j];
                    }
                }
            }
            return result;
        }

        public static Matrix<double> Subtract(this Matrix<double> first, Matrix<double> second)
        {
            int rows = first.RowCount;
            int cols = first.ColumnCount;
            var result = Matrix<double>.Build.Dense(rows, cols);

            for (int i = 0; i < rows; ++i)
            {
                for (int j = 0; j < cols; ++j)
                {
                    result[i, j] = first[i, j] - second[i, j];
                }
            }

            return result;
        }

        public static Vector<double> Sum(this Vector<double> first, Vector<double> second)
        {
            var result = Vector<double>.Build.Dense(first.Count);

            for (int i = 0; i < first.Count; ++i)
            {
                result[i] = first[i] + second[i];
            }

            return result;
        }

        public static Matrix<double> SetDiagonalElements(this Matrix<double> matrix, double value)
        {
            var n = matrix.RowCount;

            for (int i = 0; i < n; ++i)
            {
                matrix[i, i] = value;
            }

            return matrix;
        }

        public static Vector<double> GetRow(this Matrix<double> matrix, int row)
        {
            int cols = matrix.RowCount;
            var res = Vector<double>.Build.Dense(cols);

            for (int i = 0; i < cols; ++i)
            {
                res[i] = matrix[row, i];
            }

            return res;
        }

        public static Vector<double> GetColumn(this Matrix<double> matrix, int column)
        {
            int rows = matrix.RowCount;
            var res = Vector<double>.Build.Dense(rows);

            for (int i = 0; i < rows; ++i)
            {
                res[i] = matrix[i, column];
            }

            return res;
        }

        public static double SumFirstElements(this Vector<double> vector, int size)
        {
            double res = 0;
            for (int i = 0; i < size; ++i)
            {
                res += vector[i];
            }

            return res;
        }
    }
}
