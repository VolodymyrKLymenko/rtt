// Learn more about F# at http://fsharp.org

open System
open MathNet.Numerics.LinearAlgebra
open Accord.Math

let _n = 20
let _left = -3.0
let _right = 3.0
let _step = 0.0
let _function (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.00001

let mutable T = Matrix<float>.Build.Dense(_n, _n)
let mutable A = Matrix.Build.Dense(_n, _n)
let mutable B = Matrix.Build.Dense(_n, _n)
let mutable C = Matrix.Build.Dense(_n, _n)
let mutable _disceteFunction = Vector.Build.Dense(_n)

let createMatrix (v:Vector<float>) =
    let matrix = Matrix.Build.Dense(_n, _n)
    for index = 0 to (_n - 1) do
        matrix.Item(index, index) <- -2.0
        if index + 1 < _n then
            matrix.Item(index, index + 1) <- 1.0
            matrix.Item(index + 1, index) <- 1.0
    
    let vectorMatrix = Matrix.Build.Dense(_n, _n)
    for index = 0 to (_n - 1) do
        vectorMatrix.Item(index, index) <- v.Item(index)
    
    matrix * vectorMatrix

let discretizyFunction (func:(float -> float)) =
    _disceteFunction <- Vector.Build.Dense _n
    let step = (_right - _left) / ((double)_n - 1.0)
    let mutable left = _left
    for index = 0 to (_n - 1) do
        let discretizedItem = func(left)
        _disceteFunction.Item(index) <- discretizedItem
        left <- left + step
    _disceteFunction

let generateEigenValueMatrix(landa: float) =
    let vectorMatrix = Matrix.Build.Dense(_n, _n)
    for index = 0 to (_n - 1) do
        vectorMatrix.Item(index, index) <- landa
    vectorMatrix

let CreateMatrixA() =
    let result = Matrix.Build.Dense(_n, _n)
    for i = 0 to _n - 1 do
        result.[i, i] <- _disceteFunction.[i]
    let res2 = Elementwise.Multiply(T.ToArray(), result.ToArray())

    for i = 0 to _n - 1 do
        for j = 0 to _n - 1 do
            A.[i, j] <- res2.[i, j]
    A

let CreateMatrixB() =
    let result = CreateMatrixA()
    for i = 0 to _n - 1 do
        for j = 0 to _n - 1 do
            result.[i, j] <- result.[i, j] - 0.05
    B <- result
    result

let CreateMatrixC() =
    let result = Matrix.Build.Dense(_n, _n)
    for i = 0 to _n - 1 do
        for j = 0 to _n - 1 do
            result.[i, j] <- A.[i, j] - B.[i, j]
    
    C <- result
    result

[<EntryPoint>]
let main argv =
    (* let discretiziedOperator = discretizyFunction(f0)
    
    let wainshteingExecution = new WainshteinMethod(F, discretizationCount, left, right);

    wainshteingExecution.Calculate(10, out var realValues, out var calculatedValues);

    Console.WriteLine("Real values");
    foreach (var realValue in realValues.OrderByDescending(v => v))
        Console.Write($"{realValue}, ");

    Console.WriteLine("\n\n\nCalculated values");
    foreach (var calculatedValue in calculatedValues.OrderByDescending(v => v))
        Console.Write($"{calculatedValue}, ")
        *)

    Console.ReadLine();
    printfn "Hello World from F#!"
    0 // return an integer exit code
