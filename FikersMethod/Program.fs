// Learn more about F# at http://fsharp.org

open System
open MathNet.Numerics.LinearAlgebra
open Accord.Math

let discretizationVectorSize = 20
let minBorder = -1.0
let maxBorder = 1.0
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))


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
    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let res = new Decompositions.EigenvalueDecomposition(_A.ToArray(), false, true)
    let eigenValuess = res.RealEigenvalues

    Console.Write("Original values: ")
    for i = 0 to eigenValuess.Length - 1 do
        Console.Write("{0}; ", eigenValuess.[i])

    let oposite_A = _A.Inverse();

    let res = new Decompositions.EigenvalueDecomposition(oposite_A.ToArray(), false, true)
    let inversed_eigenValuess = res.RealEigenvalues

    Console.Write("\n\nCalculated values of inversed matrix: ")
    for i = 0 to eigenValuess.Length - 1 do
        Console.Write("{0}; ", inversed_eigenValuess.[i])

    let resultValues = Vector.Build.Dense(inversed_eigenValuess.Length)
    for i = 0 to inversed_eigenValuess.Length - 1 do
        resultValues.[i] <- (1.0/inversed_eigenValuess.[i])

    Console.Write("\n\nCalculated values of inversed matrix: ")
    for i = 0 to resultValues.ToArray().Length - 1 do
        Console.Write("{0}; ", resultValues.[i])

    Console.ReadKey() |> ignore
    0 // return an integer exit code
