// Learn more about F# at http://fsharp.org

open System
open MathNet.Numerics.LinearAlgebra
open Common.Extensions
open Accord.Math

let L_initial = 1.69

let discretizationVectorSize = 40
let minBorder = -1.0
let maxBorder = 1.0
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.0000001

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
        for i = r to discretizationVectorSize - 1 do

            let mutable sum = 0.0
            for j = 0 to r do
                sum <- sum + (L.[r, j] * U.[j, k])
            U.[r, k] <- operator.[r, k] - sum

            sum <- 0.0
            for j = 0 to r do
                sum <- sum + (L.[i, j] * U.[j, r])
            L.[i, r] <- (operator.[i, r] - sum) / U.[r, r]

            sum <- 0.0
            for j = 0 to r do
                sum <- sum + (M.[r, j] * U.[j, k] + L.[r, j] * V.[j, k])
            V.[r, k] <- B.[r, k] - sum

            sum <- 0.0
            for j = 0 to r do
                sum <- sum + (M.[i, j] * U.[j, r] + L.[i, j] * V.[j, r])
            M.[i, r] <- (B.[i, r] - sum - L.[i, r] * V.[r, r]) / U.[r, r]

            sum <- 0.0
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
    W <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    N <- Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)

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
    Landas.[0] <- L_initial

    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let mutable index = 1
    while index = 1 || Math.Abs(Landas.[index - 1] - Landas.[index - 2]) > percision do
        let landaMatrix = generateEigenValueMatrix(Landas.Item(index - 1))
        let D = MatrixExtensions.Subtract(_A, landaMatrix)
    
        resetMatrix()
        generateQL(D)

        let mutable sum_top = 0.0
        for k = 0 to discretizationVectorSize - 1 do
            sum_top <- sum_top + (V.[k, k] / U.[k, k])

        let mutable sum_down = 0.0
        for k = 0 to discretizationVectorSize - 1 do
            let mutable inner_sum = 0.0
            for i = 0 to discretizationVectorSize - 1 do
                if i <> k then inner_sum <- inner_sum + (V.[i, i] / U.[i, i])

            sum_down <- sum_down + (((V.[k, k] / U.[k, k]) * (V.[k, k] / U.[k, k])) - (W.[k, k] / U.[k, k])/2.0 + (V.[k, k]/U.[k, k] * inner_sum)/2.0)

        Landas.[index] <- Landas.[index - 1] - ((sum_top)/(sum_down))

        index <- index + 1

    Console.WriteLine("Helley method")
    Console.WriteLine("Amount of itterations {0}", index - 1)
    Console.WriteLine("Initial aproximation was {0}", L_initial)
    Console.WriteLine("Calculated Result eigen value: {0}\n", Landas.[index - 1])

    let inverse = new Decompositions.EigenvalueDecomposition(_A.ToArray(), false, true)
    let eigenValuess = inverse.RealEigenvalues

    Console.WriteLine ("\nReal values")
    for i = 0 to eigenValuess.Length - 1 do
        Console.Write ("{0} ,", eigenValuess.[i])

    Console.ReadKey() |> ignore
    0
