// Learn more about F# at http://fsharp.org

open System
open MathNet.Numerics.LinearAlgebra
open Common.Extensions
open Accord.Math

let L_initial = 1.9

let discretizationVectorSize = 70
let minBorder = -3.0
let maxBorder = 3.0
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.000001

let mutable U = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable V = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable L = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable M = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable N = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable W = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)


let generateQL (operator: Matrix<double>) =
    let mutable B = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    let mutable СC = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    B <- MatrixExtensions.SetDiagonalElements(B, 1.0)
    СC <- MatrixExtensions.SetDiagonalElements(СC, 1.0)

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

    for r = 0 to discretizationVectorSize - 1 do
        let mutable k = r
        for i = r + 1 to discretizationVectorSize - 1 do
            let mutable sum = 0.0
            for j = 0 to r do
                sum <- sum + (N.[r, j] * U.[j, k] + 2.0 * M.[r, j] * V.[j, k] + L.[r, j] * W.[j, k])
            W.[r, k] <- СC.[r, k] - sum

            sum <- 0.0
            for j = 0 to r do
                sum <- sum + (N.[i, j] * U.[j, r] + 2.0 * M.[i, j] * V.[j, r] + L.[i, j] * W.[j, r])
            N.[i, r] <- ((СC.[i, r] - sum  - 2.0 * M.[i, r] * V.[r, r] - L.[i, r] * W.[r, r]) / U.[r, r])

            k <- k + 1

let resetMatrix () = 
    U <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    V <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    L <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    M <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)

let createMatrix (v:Vector<float>) =
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
    let Landas = Vector.Build.Dense(1000000)
    let lower_bounds = Vector.Build.Dense(1000000)
    let higher_bounds = Vector.Build.Dense(1000000)
    let mutable calculated_landa = 0.0

    Landas.[0] <- L_initial
    lower_bounds.[0] <- minBorder
    higher_bounds.[0] <- maxBorder

    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let mutable index = 1
    while index = 1 || Math.Abs(Landas.[index - 1] - Landas.[index - 2]) > percision do
        let landaMatrix = generateEigenValueMatrix(Landas.Item(index - 1))
        let D = MatrixExtensions.Subtract(_A, landaMatrix)
        
        resetMatrix()
        generateQL(D)

        let mutable sum_min = 0.0
        let mutable sum_max = 0.0
        for k = 0 to discretizationVectorSize - 1 do
            sum_min <- sum_min + (V.[k, k] / U.[k, k])
        for k = 0 to discretizationVectorSize - 1 do
            sum_max <- sum_max + (((V.[k, k] / U.[k, k]) * (V.[k, k] / U.[k, k])) - (W.[k, k] / U.[k, k]))

        lower_bounds.[index] <- Landas.[index - 1] - abs(1.0/sum_min)
        higher_bounds.[index] <- Landas.[index - 1] - abs(sum_min/sum_max)

        Landas.[index] <- ((lower_bounds.[index] + higher_bounds.[index]) / 2.0)

        index <- index + 1

    calculated_landa <- Landas.[index - 1]

    Console.WriteLine ("\nTwo side iterations")
    Console.WriteLine ("{0,2} : {1,12}  < Result < {2,10}", "№", "m: first side", "v: second side")
    for i = 1 to index - 1 do
        Console.WriteLine ("{0,2} : {1,12}  < {2,10} < {3,10}", i, (calculated_landa - abs(lower_bounds.[i] - calculated_landa)), Landas.[i], higher_bounds.[i])

    Console.ReadKey() |> ignore
    0
