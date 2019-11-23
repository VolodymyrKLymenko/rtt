// Learn more about F# at http://fsharp.org

open MathNet.Numerics.LinearAlgebra

open Common.Extensions
open System


let pi = 3.14
let L_initial = 1.2
let q = 1.0
let N = 32
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.0001

let ii = -1.0 //Extreme.Mathematics.Complex<float>(0.0, 1.0)

let mutable U = Matrix.Build.Dense(N, N)
let mutable V = Matrix.Build.Dense(N, N)
let mutable L = Matrix.Build.Dense(N, N)
let mutable M = Matrix.Build.Dense(N, N)

let generateQL (operator: Matrix<double>) =
    let mutable B = Matrix.Build.Dense(N, N)
    B <- MatrixExtensions.SetDiagonalElements(B, 1.0)

    for r = 0 to N - 1 do
        let mutable k = r
        for i = r + 1 to N - 1 do
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
    U <- Matrix.Build.Dense(N, N)
    V <- Matrix.Build.Dense(N, N)
    L <- Matrix.Build.Dense(N, N)
    M <- Matrix.Build.Dense(N, N)

let createMatrix (v:Vector<float>) =
    let matrix = Matrix.Build.Dense(N, N)
    for index = 0 to (N - 1) do
        matrix.Item(index, index) <- -2.0
        if index + 1 < N then
            matrix.Item(index, index + 1) <- 1.0
            matrix.Item(index + 1, index) <- 1.0
    
    let vectorMatrix = Matrix.Build.Dense(N, N)
    for index = 0 to (N - 1) do
        vectorMatrix.Item(index, index) <- v.Item(index)
    
    matrix * vectorMatrix

let discretizyFunction (func:(float -> float)) =
    let resVector = Vector.Build.Dense N
    let step = (-q - q) / ((double)N - 1.0)
    let mutable left = -q
    for index = 0 to (N - 1) do
        let discretizedItem = func(left)
        resVector.Item(index) <- discretizedItem
        left <- left + step
    resVector

let generateEigenValueMatrix(landa: float) =
    let vectorMatrix = Matrix.Build.Dense(N, N)
    for index = 0 to (N - 1) do
        vectorMatrix.Item(index, index) <- landa
    vectorMatrix

let pow(x: double, n: int) =
    let mutable res = 1.0
    for i = 0 to n - 1 do
        res <- res * x
    res

let landa_j(i: double, j: double) =
    (L_initial + q + exp((double)ii * 2.0 * j / (double)N))

[<EntryPoint>]
let main argv =
    let mutable m = 0.0

    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let landaMatrix = generateEigenValueMatrix(L_initial)
    let D = MatrixExtensions.Subtract(_A, landaMatrix)
    generateQL(D)

    let mutable del = 0.0
    for r = 0 to N - 1 do
        del <- del + (V.[r, r] / U.[r, r])

    let m = abs((int)((double)(1.0 / (2.0 * pi * ii)) * del))

    let landas = Array.create m  0.0

    let s = Vector<double>.Build.Dense(m)

    for k = 0 to m - 1 do
        let mutable sum = 0.0
        for j = 0 to N - 1 do
            let landa = L_initial + (q * exp((double)k * (2.0 * pi * (double)j) / (double)N))
            let landaMatrix = generateEigenValueMatrix(landa)
            let D = MatrixExtensions.Subtract(_A, landaMatrix)
            generateQL(D)

            let mutable del = 0.0
            for r = 0 to N - 1 do
                del <- del + (V.[r, r] / U.[r, r])

            sum <- sum + (pow(landa, k) * q * exp((double)ii * (2.0 * pi * (double)j / (double)N)) * del)
        
        s.[k] <- (1.0 / (double)m) * sum

    for j = 1 to m do
        landas.[j - 1] <- L_initial + q * exp((double(ii) * 2.0 * pi * (double)j / (double)m))

    for j = 0 to (m - 1) do
        let mutable previousValue = landas.[j]
        let mutable index = 1
        while (index = 1 || abs(previousValue - landas.[j]) > percision) do
            let F = Vector<double>.Build.Dense(m)

            for k = 0 to (m - 1) do
                let mutable phi_k = 0.0
                for p = 0 to m - 1 do
                    phi_k <- phi_k + pow(landa_j((double)k, (double)p), k)
                F.[k] <- phi_k
                   
            previousValue <- landas.[j]
            landas.[j] <- landas.[j] - (1.0 / (MatrixExtensions.MultiplyVectorsTansponential(F, MatrixExtensions.SubtractVectors(F, s))))
            index <- index + 1

            if index < 100 then
                landas.[j] <- -1.0

    for i = 0 to (m - 1) do
        if landas.[i] > 0.0 then Console.WriteLine ("L_{0} :: {1}", i + 1, landas.[i])

    printfn "End of execution"
    Console.ReadKey() |> ignore
    0 // return an integer exit code
