// Learn more about F# at http://fsharp.org

open System

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra.Double

open Extreme.Mathematics

open  Common.Extensions

let discretizationVectorSize = 10
let minBorder = -1.0
let maxBorder = 1.0
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.0001

let mutable U = LinearAlgebra.Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable V = LinearAlgebra.Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable L = LinearAlgebra.Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable M = LinearAlgebra.Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)


let generateQL (operator: LinearAlgebra.Matrix<double>) =
    let mutable B = LinearAlgebra.Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    B <- MatrixExtensions.SetDiagonalElements(B, 1.0)

    for r = 0 to discretizationVectorSize - 1 do
        let mutable k = r
        for i = r + 1 to discretizationVectorSize - 1 do
            let mutable L_row = MatrixExtensions.GetRow(L, k)
            let mutable U_column = MatrixExtensions.GetColumn(U, k)
            U.Item(r, k) <- operator.Item(r, k) - MatrixExtensions.SumFirstElements(MatrixExtensions.MultiplyVectors(L_row, U_column), r)

            L_row <- MatrixExtensions.GetRow(L, i)
            U_column <- MatrixExtensions.GetColumn(U, r)
            L.Item(i, r) <- (operator.Item(i, r) - MatrixExtensions.SumFirstElements(MatrixExtensions.MultiplyVectors(L_row, U_column), r)) / U.Item(r, r)

            let mutable sum_1 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(M, r), MatrixExtensions.GetColumn(U, k));
            let mutable sum_2 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, r), MatrixExtensions.GetColumn(V, k))
            V.Item(r, k) <- B.Item(r, k) - MatrixExtensions.SumFirstElements(MatrixExtensions.Sum(sum_1, sum_2), r)
            
            sum_1 <- MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(M, i), MatrixExtensions.GetColumn(U, r))
            sum_2 <- MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, i), MatrixExtensions.GetColumn(V, r))
            M.Item(i, r) <- (B.Item(i, r) - MatrixExtensions.SumFirstElements(MatrixExtensions.Sum(sum_1, sum_2), r) - L.[i, r] * V.[r, r]) / U.[r, r]

            k <- k + 1

        k <- discretizationVectorSize - 1

        U.[r, k] <-  operator.[r, k] - MatrixExtensions.SumFirstElements(MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, r), MatrixExtensions.GetColumn(U, k)), r)
        let sum_1 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(M, r), MatrixExtensions.GetColumn(U, k))
        let sum_2 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, r), MatrixExtensions.GetColumn(V, k))
        V.[r, k] <-  B.[r, k] - MatrixExtensions.SumFirstElements(MatrixExtensions.Sum(sum_1, sum_2), r)

let resetMatrix () = 
    U <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    V <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    L <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    M <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)

let createMatrix (v:LinearAlgebra.Vector<float>) =
    let matrix = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    for index = 0 to (discretizationVectorSize - 1) do
        matrix.Item(index, index) <- -2.0
        if index + 1 < discretizationVectorSize then
            matrix.Item(index, index + 1) <- 1.0
            matrix.Item(index + 1, index) <- 1.0
    
    let vectorMatrix = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    for index = 0 to (discretizationVectorSize - 1) do
        vectorMatrix.Item(index, index) <- v.Item(index)
    
    matrix * vectorMatrix

let discretizyFunction (func:(float -> float)) =
    let resVector = Vector.Build.Dense discretizationVectorSize
    let step = (maxBorder - minBorder) / ((double)discretizationVectorSize - 1.0)
    let mutable left = minBorder
    for index = 0 to (discretizationVectorSize - 1) do
        let discretizedItem = func(left)
        resVector.Item(index) <- discretizedItem
        left <- left + step
    resVector

let generateEigenValueMatrix(landa: float) =
    let vectorMatrix = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    for index = 0 to (discretizationVectorSize - 1) do
        vectorMatrix.Item(index, index) <- landa
    vectorMatrix

[<EntryPoint>]
let main argv =
    let dL_values = Vector.Build.Dense(10000)
    let eigenvalues = Vector.Build.Dense(10000)
    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    eigenvalues.Item(0) <- 0.4
    Console.WriteLine("Iteration - {0,2}   :: Next landa value is " + (eigenvalues.Item(0)).ToString(), "1")

    let mutable index = 1
    while index = 1 || Math.Abs(dL_values.Item(index) - dL_values.Item(index - 1)) > percision do
        let landaMatrix = generateEigenValueMatrix(eigenvalues.Item(index - 1))
        let D = MatrixExtensions.Subtract(_A, landaMatrix)
        
        generateQL(D)

        let mutable sum = 0.0
        for i = 0 to discretizationVectorSize - 1 do
            sum <- sum + (V.Item(i, i) / U.Item(i, i))
        dL_values.Item(index) <- (1.0 / sum)

        eigenvalues.Item(index) <- (eigenvalues.Item(index - 1) + dL_values.Item(index))

        Console.WriteLine("Iteration - {0,2}   :: Next landa value is " + (eigenvalues.Item(index)).ToString(), (index + 1).ToString())

        index <- index + 1
        resetMatrix()

    Console.ReadKey() |> ignore

    printfn "Hello World from F#!"
    0 // return an integer exit code
