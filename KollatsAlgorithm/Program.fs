open System

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra.Double

let mutable discretizationVectorSize = 10
let minBorder = -1.0
let maxBorder = 1.0
let l2 = 0.7
let e = 0.000001

let executeScalarProduct (v1: LinearAlgebra.Vector<float>, v2: LinearAlgebra.Vector<float>) =
    let mutable res = 0.0
    for i = 0 to (discretizationVectorSize - 1) do
        res <- res + (v1.[i] * v2.[i])
    res

let discretizyFunction (func:(float -> float)) =
    let resVector = Vector.Build.Dense discretizationVectorSize
    let step = (maxBorder - minBorder) / ((double)discretizationVectorSize - 1.0)
    let mutable left = minBorder
    for index = 0 to (discretizationVectorSize - 1) do
        let discretizedItem = func(left)
        resVector.Item(index) <- discretizedItem
        left <- left + step
    resVector

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
    
    vectorMatrix * matrix

[<EntryPoint>]
let main argv =
    // Initialization
    discretizationVectorSize <-10
    let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))

    // Preparation
    let f0_discr = discretizyFunction(f0)

    let T = createMatrix(f0_discr)
    let a0 = executeScalarProduct(f0_discr, f0_discr)

    // Iterations
    let m_valus = Array.create 10000 0.0
    let v_valus = Array.create 10000 0.0
    let landa_valus = Array.create 10000 0.0
    
    let mutable f_k_minus_1_discretizied = f0_discr
    let mutable index = 1

    Console.WriteLine ("\nKollats algorithm iterations")
    Console.WriteLine ("{0,3} {1,10} {2,15} {3,15}", "Iter.", 'v', "landa", "m")

    while index = 1 || Math.Abs(landa_valus.[index - 1] - landa_valus.[index - 2]) > e do
        let fk_discr = T.Solve(f_k_minus_1_discretizied)
        let a_k = executeScalarProduct(fk_discr, f0_discr)
        
        if index = 1 then
            let m_0 = a_k / a0
            Array.set m_valus 0 m_0

        f_k_minus_1_discretizied <- fk_discr

        let fk_plus_1_discr = T.Solve(fk_discr)
        let a_k_plus_1 = executeScalarProduct(fk_plus_1_discr, f0_discr)

        let m_k = m_valus.[index - 1]
        let m_k_plus_1 = a_k / a_k_plus_1

        Array.set m_valus index m_k_plus_1

        Console.WriteLine ("{0,2}|\t{1:F10}\t{2:F10}\t{3:F10}", index, v_valus.[index - 1], landa_valus.[index - 1], m_valus.[index - 1])

        let v_k_plus_1 = m_k_plus_1 - (((m_k - m_k_plus_1) * m_k_plus_1) / (l2 - m_k_plus_1))
        let landa_k_plus_1 = (m_k_plus_1 + v_k_plus_1) / 2.0

        Array.set v_valus index v_k_plus_1
        Array.set landa_valus index landa_k_plus_1

        index <- index + 1

    Console.ReadKey() |> ignore

    Console.WriteLine ("\nKrilov borders")
    Console.WriteLine ("{0,12}  < Result < {1,10}", "    m lwer    ", "m upper")

    for index = 1 to (index - 1) do
        let krilov_min = m_valus.[index] - Math.Sqrt(Math.Abs(m_valus.[index + 1] - m_valus.[index]) * m_valus.[index + 1])
        let krilov_max = m_valus.[index] + Math.Sqrt(Math.Abs(m_valus.[index + 1] - m_valus.[index]) * m_valus.[index + 1])

        Console.WriteLine ("{0,17} < L < {1,17}", krilov_min.ToString(), krilov_max.ToString())

    Console.ReadKey() |> ignore
    0 // return an integer exit code