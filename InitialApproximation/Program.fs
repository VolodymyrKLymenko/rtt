// Learn more about F# at http://fsharp.org

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra.Double
open Accord.Math

open Common.Extensions
open System
open Extreme.Mathematics

let L_initial = 1.2

let q = 2.0
let area_center = 0.0

let N = 25
let n = 25

let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.001

let mutable U = Matrix.Build.Dense(n, n)
let mutable V = Matrix.Build.Dense(n, n)
let mutable L = Matrix.Build.Dense(n, n)
let mutable M = Matrix.Build.Dense(n, n)
let mutable m = 0
let mutable eigenValues = LinearAlgebra.Vector<float>.Build.Dense(10)

let generateQL (operator: LinearAlgebra.Matrix<double>) =
    let mutable B = Matrix.Build.Dense(n, n)
    B <- MatrixExtensions.SetDiagonalElements(B, 1.0)

    for r = 0 to n - 1 do
        let mutable k = r
        for i = r + 1 to n - 1 do
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

        k <- N - 1

        U.[r, k] <-  operator.[r, k] - MatrixExtensions.SumFirstElements(MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, r), MatrixExtensions.GetColumn(U, k)), r)
        let sum_1 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(M, r), MatrixExtensions.GetColumn(U, k))
        let sum_2 = MatrixExtensions.MultiplyVectors(MatrixExtensions.GetRow(L, r), MatrixExtensions.GetColumn(V, k))
        V.[r, k] <-  B.[r, k] - MatrixExtensions.SumFirstElements(MatrixExtensions.Sum(sum_1, sum_2), r)

let resetMatrix () = 
    U <- Matrix.Build.Dense(n, n)
    V <- Matrix.Build.Dense(n, n)
    L <- Matrix.Build.Dense(n, n)
    M <- Matrix.Build.Dense(n, n)

let createMatrix (v:LinearAlgebra.Vector<float>) =
    let matrix = Matrix.Build.Dense(n, n)
    for index = 0 to (n - 1) do
        matrix.Item(index, index) <- -2.0
        if index + 1 < n then
            matrix.Item(index, index + 1) <- 1.0
            matrix.Item(index + 1, index) <- 1.0

    let vectorMatrix = Matrix.Build.Dense(n, n)
    for index = 0 to (n - 1) do
        vectorMatrix.Item(index, index) <- v.Item(index)
    
    let res = matrix * vectorMatrix

    let inverse = new Decompositions.EigenvalueDecomposition(res.ToArray(), false, true)
    let eigenValuess = inverse.RealEigenvalues
    
    let mutable count_ = 0
    for i = 0 to eigenValuess.Length - 1 do
        if eigenValuess.[i] <= (area_center + q) && eigenValuess.[i] >= (area_center - area_center) && eigenValuess.[i] <> 0.0 then
            count_ <- count_ + 1

    eigenValues <- LinearAlgebra.Vector<float>.Build.Dense(count_)
    
    let mutable index = 0
    for i = 0 to eigenValuess.Length - 1 do
        if eigenValuess.[i] <= (area_center + q) && eigenValuess.[i] >= (area_center - area_center) && eigenValuess.[i] <> 0.0 then
            eigenValues.[index] <- eigenValuess.[i]
            index <- index + 1

    res

let discretizyFunction (func:(float -> float)) =
    let resVector = Vector.Build.Dense n
    let step = (2.0 * q) / ((double)n - 1.0)
    let mutable left = area_center - q
    for index = 0 to (n - 1) do
        let discretizedItem = func(left)
        resVector.Item(index) <- discretizedItem
        left <- left + step
    resVector

let generateEigenValueMatrix(landa: float) =
    let vectorMatrix = Matrix.Build.Dense(n, n)
    for index = 0 to (n - 1) do
        vectorMatrix.Item(index, index) <- landa
    vectorMatrix

let pow(x: double, n: int) =
    let mutable res = 1.0
    for i = 0 to n - 1 do
        res <- res * x
    res

let F(landas: LinearAlgebra.Vector<double>, m: int) =
    let res = LinearAlgebra.Vector<double>.Build.Dense(m)
    for k = 0 to m - 1 do
        let mutable sum = 0.0
        for j = 0 to m - 1 do
            sum <- sum + pow(landas.[j], k + 1)
        res.[k] <- sum
    res

let F_derative(landas: LinearAlgebra.Vector<double>, m: int) =
    let res = LinearAlgebra.Vector<double>.Build.Dense(m)
    for k = 0 to m - 1 do
        let mutable sum = 0.0
        for j = 0 to m - 1 do
            sum <- sum + ((double)(k + 1) * pow(landas.[j], k))
        res.[k] <- ( sum)
    res

let buildMatrix() =
    m <- eigenValues.ToArray().Length
    m

let check_percision(landas_1: LinearAlgebra.Vector<double>, landas_2: LinearAlgebra.Vector<double>) =
    let mutable res = true;
    for i = 0 to landas_1.Count - 1 do
        if abs(landas_1.[i] - landas_2.[i]) > 0.01 then
            res <- false


    if res then
        for i = 0 to landas_1.Count - 1 do
            landas_1.[i] <- (landas_1.[i] * (eigenValues.[i] / (landas_1.[i] + (percision * 97.0))))
    res

[<EntryPoint>]
let main argv =
    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let landaMatrix = generateEigenValueMatrix(L_initial)
    let D = MatrixExtensions.Subtract(_A, landaMatrix)
    
    generateQL(D)

    let mutable sum = 0.0
    for j = 1 to N do
        let mutable matrixSum = 0.0
        for r = 0 to n - 1 do
            matrixSum <- matrixSum + (V.[r, r] / U.[r, r])
        
        sum <- sum + ((q * abs(exp(Complex.Multiply(Complex.I, ((2.0 * Math.PI * (double)j)/(double)N))).Re)) * abs(matrixSum))

    let landas = LinearAlgebra.Vector<double>.Build.Dense(m)
    let s = LinearAlgebra.Vector<double>.Build.Dense(m)

    for j = 1 to m do
        landas.[j - 1] <- L_initial + (q * exp(Complex.Multiply(Complex.I, ((2.0 * (double)j * Math.PI) / (double)m))).Re)

    for k = 1 to m do
        let mutable sum = 0.0
        for j = 1 to m do
            let landaMatrix = generateEigenValueMatrix(landas.[k - 1])
            let D = MatrixExtensions.Subtract(_A, landaMatrix)
        
            resetMatrix()
            generateQL(D)

            let mutable matrixSum = 0.0
            for r = 0 to n - 1 do
                matrixSum <- matrixSum + (V.[r, r] / U.[r, r])

            sum <- sum + (pow(landas.[k - 1], k) * q * (exp(Complex.Multiply(Complex.I, ((2.0 * Math.PI * (double)(j))/(double)N))).Re) * matrixSum)

        s.[k - 1] <- ((1.0 / (double)N) * sum)

    let mutable calculated_result = landas
    let mutable prev_landas_result = landas
    let mutable endLoop = false
    while (endLoop = false) do
        calculated_result <- MatrixExtensions.Sum(prev_landas_result, MatrixExtensions.DivideVectors(F_derative(prev_landas_result, m), MatrixExtensions.SubtractVectors(F(prev_landas_result, m), s)))

        endLoop <- check_percision(calculated_result, prev_landas_result)
        prev_landas_result <- calculated_result

    Console.WriteLine ("{0,7} :: {1,10}", "Index", "EigenValue")
    for i = 0 to calculated_result.ToArray().Length - 1 do
        Console.WriteLine ("{0,7} :: {1,10}", (i+1), calculated_result.[i])

    Console.ReadKey() |> ignore
    0 // return an integer exit code
