namespace Common.Extensions

open MathNet.Numerics

open System.Runtime.CompilerServices

module MatrixExtensions =
    let MultiplyVectors(first: LinearAlgebra.Vector<double>, second: LinearAlgebra.Vector<double>) =
        let result = LinearAlgebra.Vector<double>.Build.Dense(first.Count)

        for i = 0 to first.Count - 1 do
            result.Item(i) <- first.[i] * second.[i]

        result

    let DivideVectors(first: LinearAlgebra.Vector<double>, second: LinearAlgebra.Vector<double>) =
        let result = LinearAlgebra.Vector<double>.Build.Dense(first.Count)

        for i = 0 to first.Count - 1 do
            result.Item(i) <- first.[i] / second.[i]

        result

    let MultiplyVectorsTansponential(first: LinearAlgebra.Vector<double>, second: LinearAlgebra.Vector<double>) =
        let mutable result = 0.0

        for i = 0 to first.Count - 1 do
            result <- result + (first.[i] * second.[i])

        result

    let MultiplyMatrix(first: LinearAlgebra.Matrix<double>, second: LinearAlgebra.Matrix<double>) =
        let rows = first.RowCount
        let cols = first.ColumnCount
        let result = LinearAlgebra.Matrix<double>.Build.Dense(rows, cols)
        
        for i = 0 to rows - 1 do
            for j = 0 to cols - 1 do
                for k = 0 to cols - 1 do
                    result.Item(i, j) <- result.[i, j] + first.[i, k] * second.[k, j]

        result

    let Subtract(first: LinearAlgebra.Matrix<double>, second: LinearAlgebra.Matrix<double>) =
        let rows = first.RowCount
        let cols = first.ColumnCount
        let result = LinearAlgebra.Matrix<double>.Build.Dense(rows, cols)
        
        for i = 0 to rows - 1 do
            for j = 0 to cols - 1 do
                result.Item(i, j) <- first.[i, j] - second.[i, j]
        
        result

    let SubtractVectors(first: LinearAlgebra.Vector<double>, second: LinearAlgebra.Vector<double>) =
        let rows = first.Count
        let result = LinearAlgebra.Vector<double>.Build.Dense(rows)
    
        for i = 0 to rows - 1 do
            result.Item(i) <- first.[i] - second.[i]
    
        result

    let SubtractValue(first: LinearAlgebra.Vector<double>, value: double) =
        let rows = first.Count
        let result = LinearAlgebra.Vector<double>.Build.Dense(rows)

        for i = 0 to rows - 1 do
            result.Item(i) <- first.[i] - value

        result

    let InverseVector(vector: LinearAlgebra.Vector<double>) =
        let rows = vector.Count
        let result = LinearAlgebra.Vector<double>.Build.Dense(rows)

        for i = 0 to rows - 1 do
            result.Item(i) <- 1.0 / vector.[i]

        result

    let Sum(first: LinearAlgebra.Vector<double>, second: LinearAlgebra.Vector<double>) =
        let result = LinearAlgebra.Vector<double>.Build.Dense(first.Count)
        
        for i = 0 to first.Count - 1 do
            result.[i] <- first.[i] + second.[i]
        
        result

    let SetDiagonalElements(matrix: LinearAlgebra.Matrix<double>, value: double) =
        let n = matrix.RowCount

        for i = 0 to n - 1 do
            matrix.[i, i] <- value

        matrix

    let GetRow(matrix: LinearAlgebra.Matrix<double>, row: int) =
        let cols = matrix.RowCount
        let res = LinearAlgebra.Vector<double>.Build.Dense(cols)

        for i = 0 to cols - 1 do
            res.[i] <- matrix.[row, i]

        res

    let GetColumn(matrix: LinearAlgebra.Matrix<double>, column: int) =
        let rows = matrix.RowCount
        let res = LinearAlgebra.Vector<double>.Build.Dense(rows)

        for i = 0 to rows - 1 do
            res.[i] <- matrix.[i, column]

        res

    let SumFirstElements(vector: LinearAlgebra.Vector<double>, size: int) =
        let mutable res = 0.0;

        for i = 0 to size - 1 do
            res <- res + vector.[i]

        res

[<Extension>]
type MatrixExtensionsss () =
    [<Extension>]
    static member SumFirstElements(ref: Ref<LinearAlgebra.Vector<double>>, size: int) =
        let mutable res = 0.0;
    
        for i = 0 to size do
            res <- res + ref.Value.[i]

        res

    [<Extension>]
    static member MultiplyVectors(first: Ref<LinearAlgebra.Vector<double>>, second: LinearAlgebra.Vector<double>) =
        let result = LinearAlgebra.Vector<double>.Build.Dense(first.Value.Count)
        
        for i = 0 to first.Value.Count do
            result.Item(i) <- first.Value.[i] * second.[i]
        
        result
    
    [<Extension>]
    static member MultiplyMatrix(first: Ref<LinearAlgebra.Matrix<double>>, second: LinearAlgebra.Matrix<double>) =
        let rows = first.Value.RowCount
        let cols = first.Value.ColumnCount
        let result = LinearAlgebra.Matrix<double>.Build.Dense(rows, cols)

        for i = 0 to rows do
            for j = 0 to cols do
                for k = 0 to cols do
                    result.Item(i, j) <- result.[i, j] + first.Value.[i, k] * second.[k, j]

        result

    [<Extension>]
    static member Subtract(first: Ref<LinearAlgebra.Matrix<double>>, second: LinearAlgebra.Matrix<double>) =
        let rows = first.Value.RowCount
        let cols = first.Value.ColumnCount
        let result = LinearAlgebra.Matrix<double>.Build.Dense(rows, cols)

        for i = 0 to rows do
            for j = 0 to cols do
                result.Item(i, j) <- first.Value.[i, j] - second.[i, j]

        result

    [<Extension>]
    static member Sum(first: Ref<LinearAlgebra.Vector<double>>, second: LinearAlgebra.Vector<double>) =
        let result = LinearAlgebra.Vector<double>.Build.Dense(first.Value.Count)

        for i = 0 to first.Value.Count do
            result.[i] <- first.Value.[i] + second.[i]

        result

    [<Extension>]
    static member SetDiagonalElements(matrix: Ref<LinearAlgebra.Matrix<double>>, value: double) =
        let n = matrix.Value.RowCount
        
        for i = 0 to n do
            matrix.Value.[i, i] <- value
        
        matrix.Value

    [<Extension>]
    static member GetRow(matrix: Ref<LinearAlgebra.Matrix<double>>, row: int) =
        let cols = matrix.Value.RowCount
        let res = LinearAlgebra.Vector<double>.Build.Dense(cols)
        
        for i = 0 to cols do
            res.[i] <- matrix.Value.[row, i]
        
        res