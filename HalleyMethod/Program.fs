// Learn more about F# at http://fsharp.org

open System
open MathNet.Numerics.LinearAlgebra
open Common.Extensions
open Accord.Math

let L_initial = 0.01

let discretizationVectorSize = 20
let minBorder = -3.0
let maxBorder = 3.0
let f0 (x) = -((4.0 - x * x) * (4.0 - x * x))
let percision = 0.00001

let mutable U = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable V = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable L = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable M = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable N = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
let mutable W = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)

let generateQL (operator: Matrix<double>) =
    let mutable B = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    let mutable CC = Matrix.Build.Dense(discretizationVectorSize, discretizationVectorSize)
    B <- MatrixExtensions.SetDiagonalElements(B, 1.0)
    CC <- MatrixExtensions.SetDiagonalElements(CC, 1.0)

    for r = 0 to discretizationVectorSize - 1 do
        let mutable k = r
        for i = r + 1 to discretizationVectorSize - 1 do
            let mutable sum = 0.0
            for j = 0 to r - 1 do
                sum <- sum + (L.[r, j] * U.[j, k])
            U.[r, k] <- operator.[r, k] - sum

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (L.[i, j] * U.[j, r])
            L.[i, r] <- (operator.[i, r] - sum) / U.[r, r]

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (M.[r, j] * U.[j, k] + L.[r, j] * V.[j, k])
            V.[r, k] <- B.[r, k] - sum

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (M.[i, j] * U.[j, r] + L.[i, j] * V.[j, r])
            M.[i, r] <- (B.[i, r] - sum - L.[i, r] * V.[r, r]) / U.[r, r]

            B <-  M * U + L * V

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (N.[r, j] * U.[j, k] + 2.0 * M.[r, j] * V.[j, k] + L.[r, j] * W.[j, k])
            W.[r, k] <- CC.[r, k] - sum

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (N.[i, j] * U.[j, r] + 2.0 * M.[i, j] * V.[j, r] + L.[i, j] * W.[j, r])
            N.[i, r] <- ((CC.[i, r] - sum  - 2.0 * M.[i, r] * V.[r, r] - L.[i, r] * W.[r, r]) / U.[r, r])

            k <- k + 1

        if k <= discretizationVectorSize - 1 then
            let mutable sum = 0.0
            for j = 0 to r - 1 do
                sum <- sum + (L.[r, j] * U.[j, k])
            U.[r, k] <- operator.[r, k] - sum

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (M.[r, j] * U.[j, k] + L.[r, j] * V.[j, k])
            V.[r, k] <- B.[r, k] - sum

            B <-  M * U + L * V

            sum <- 0.0
            for j = 0 to r - 1 do
                sum <- sum + (N.[r, j] * U.[j, k] + 2.0 * M.[r, j] * V.[j, k] + L.[r, j] * W.[j, k])
            W.[r, k] <- CC.[r, k] - sum

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

let pow(x: double, n: int) =
    let mutable res = 1.0
    for i = 0 to n - 1 do
        res <- res * x
    res

let F_landa() =
    let mutable result = 1.0
    for i = 0 to discretizationVectorSize - 1 do
        result <- result * U.[i, i]
    result

let F_second_derative_landa() =
    let mutable left_sum = 0.0
    let mutable right_sum = 0.0

    for k = 0 to discretizationVectorSize - 1 do
        let mutable u_exclude_k_prod = 1.0
        for i = 0 to discretizationVectorSize - 1 do
            if i <> k then u_exclude_k_prod <- u_exclude_k_prod * U.[i, i]
        left_sum <- left_sum + (W.[k, k] * u_exclude_k_prod)

    for k = 0 to discretizationVectorSize - 1 do
        let mutable v_exclude_k_sum = 0.0
        for j = 0 to discretizationVectorSize - 1 do
            let mutable u_exclude_k_j_prod = 1.0
            for i = 0 to discretizationVectorSize - 1 do
                if i <> k && i <> j then u_exclude_k_j_prod <- u_exclude_k_j_prod * U.[i, i]
            
            if j <> k then v_exclude_k_sum <- v_exclude_k_sum + (V.[j, j] * u_exclude_k_j_prod)
        right_sum <- right_sum + (V.[k, k] * v_exclude_k_sum)

    left_sum + right_sum

[<EntryPoint>]
let main argv =
    let mutable Landas = Vector.Build.Dense(1000000)
    Landas.[0] <- L_initial
    let mutable result = 0.0

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

    result <- Landas.[index - 1]
    Console.WriteLine("One side Helley method (3rd count)")
    Console.WriteLine("Initial aproximation was {0}", L_initial)
    Console.WriteLine("Calculated Result eigen value: {0}\n", result)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Landas <- Vector.Build.Dense(100000)
    let mutable M_values = Vector<double>.Build.Dense(100000)
    Landas.[0] <- L_initial
    M_values.[0] <- L_initial

    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    let f_L = F_landa()
    let f_second_derative_L = F_second_derative_landa()

    index <- 1
    while index = 1 || Math.Abs(Landas.[index - 1] - M_values.[index - 1]) > percision do
        let landa = Landas.[index - 1]
        let landaMatrix = generateEigenValueMatrix(landa)
        let D = MatrixExtensions.Subtract(_A, landaMatrix)
    
        resetMatrix()
        generateQL(D)

        if f_L * f_second_derative_L < 0.0 then
            let mutable f_m_landa = 0.0
            let mutable f_m_derative_landa = 0.0
            let mutable h_m_derative_landa = 0.0

            for k = 0 to discretizationVectorSize - 1 do
                f_m_landa <- f_m_landa + (V.[k, k] / U.[k, k])
                f_m_derative_landa <- f_m_derative_landa + (((V.[k, k]/U.[k, k])*(V.[k, k]/U.[k, k])) - (W.[k, k]/U.[k, k]))

                let mutable inner_sum = 0.0
                for i = 0 to discretizationVectorSize - 1 do
                    if i <> k then inner_sum <- inner_sum + (V.[i, i]/U.[i, i])
                h_m_derative_landa <- h_m_derative_landa + (((V.[k, k]/U.[k, k])*(V.[k, k]/U.[k, k])) - ((1.0/2.0)*(W.[k, k]/U.[k, k])) + ((1.0/2.0)*(V.[k, k]/U.[k, k])*inner_sum))

            M_values.[index] <- Landas.[index - 1] - ((f_m_landa)/(f_m_derative_landa))
            Landas.[index] <- Landas.[index - 1] - ((f_m_landa)/(h_m_derative_landa))
        else
            let mutable f_m_landa = 0.0
            let mutable f_m_derative_landa = 0.0
            let mutable h_m_derative_landa = 0.0

            for k = 0 to discretizationVectorSize - 1 do
                f_m_landa <- f_m_landa + ((V.[k, k] / U.[k, k]))
                f_m_derative_landa <- f_m_derative_landa + (((V.[k, k]/U.[k, k])*(V.[k, k]/U.[k, k])) - W.[k, k]/U.[k, k])

                let mutable inner_sum = 0.0
                for i = 0 to discretizationVectorSize - 1 do
                    if i <> k then inner_sum <- inner_sum + (V.[i, i]/U.[i, i])
                h_m_derative_landa <- h_m_derative_landa + (((V.[k, k]/U.[k, k])*(V.[k, k]/U.[k, k])) - (1.0/2.0)*(W.[k, k]/U.[k, k]) + (1.0/2.0)*(V.[k, k]/U.[k, k])*inner_sum)

            M_values.[index] <- Landas.[index - 1] - (1.0/f_m_derative_landa)
            Landas.[index] <- Landas.[index - 1] - (f_m_landa/h_m_derative_landa)

        index <- index + 1

    Console.WriteLine("\nTwo side inclusive Helley method 1")
    Console.WriteLine("Initial aproximation was {0}", L_initial)
    Console.WriteLine("Calculated Result eigen value:")
    let res = Landas.[index - 1];
    let total = index - 1
    let step = total / 14
    index <- 0
    for i = 0 to total do
        if i = index || i = total then
            Console.WriteLine("     -> Landa: {0} < L < M: {1}", Landas.[i] - (abs(res - Landas.[i])), M_values.[i] + (abs(res - Landas.[i])))
            index <- index + step

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    let v1 = Landas.[1] - (abs(res - Landas.[1]))
    let m2 = M_values.[1] + (abs(res - Landas.[1]))
    Landas <- Vector.Build.Dense(100000)
    M_values <- Vector<double>.Build.Dense(100000)
    Landas.[0] <- v1
    M_values.[0] <- m2

    let discretiziedOperator = discretizyFunction(f0)
    let _A = createMatrix(discretiziedOperator)

    index <- 1
    while index = 1 || Math.Abs(Landas.[index - 1] - M_values.[index - 1]) > percision do
        let mutable landaMatrix = generateEigenValueMatrix(Landas.[index - 1])
        let mutable D = MatrixExtensions.Subtract(_A, landaMatrix)
    
        resetMatrix()
        generateQL(D)

        let mutable f_m_landa = 0.0
        let mutable h_m_derative_landa = 0.0

        for k = 0 to discretizationVectorSize - 1 do
            f_m_landa <- f_m_landa + (V.[k, k] / U.[k, k])

            let mutable inner_sum = 0.0
            for i = 0 to discretizationVectorSize - 1 do
                if i <> k then inner_sum <- inner_sum + (V.[i, i]/U.[i, i])
            h_m_derative_landa <- h_m_derative_landa + (pow((V.[k, k]/U.[k, k]), 2) - (1.0/2.0)*(W.[k, k]/U.[k, k]) + (1.0/2.0)*(V.[k, k]/U.[k, k])*inner_sum)

        Landas.[index] <- Landas.[index - 1] - (f_m_landa/h_m_derative_landa)


        landaMatrix <- generateEigenValueMatrix(M_values.[index - 1])
        D <- MatrixExtensions.Subtract(_A, landaMatrix)
    
        resetMatrix()
        generateQL(D)

        let mutable f_m_landa = 0.0
        let mutable h_m_derative_landa = 0.0

        for k = 0 to discretizationVectorSize - 1 do
            f_m_landa <- f_m_landa + (V.[k, k] / U.[k, k])

            let mutable inner_sum = 0.0
            for i = 0 to discretizationVectorSize - 1 do
                if i <> k then inner_sum <- inner_sum + (V.[i, i]/U.[i, i])
            h_m_derative_landa <- h_m_derative_landa + (pow((V.[k, k]/U.[k, k]), 2) - (1.0/2.0)*(W.[k, k]/U.[k, k]) + (1.0/2.0)*(V.[k, k]/U.[k, k])*inner_sum)

        M_values.[index] <- Landas.[index - 1] - (f_m_landa/h_m_derative_landa)

        index <- index + 1

    Console.WriteLine("\nTwo side inclusive Helley method 2")
    Console.WriteLine("Initial aproximation was {0}", L_initial)
    for i = 1 to index - 1 do
        Console.WriteLine("     -> Landa: {0} < L < M: {1}", Landas.[i], M_values.[i])

    Console.ReadKey() |> ignore
    0
