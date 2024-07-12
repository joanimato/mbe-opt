namespace Optimize

// #r "nuget: FSharp.Stats"

// #r "nuget: MathNet.Numerics, 5.0.0-alpha08"
// #r "nuget: MathNet.Numerics.FSharp, 5.0.0-alpha08"



module GradientDescent = 

    open Output.ReadOutput
    open Input.InputMBE
    // open FSharp.Stats 
    // open Extreme.Mathematics.Optimization
    // open MathNet.Numerics.Optimization 
    open MathNet.Numerics.LinearAlgebra

    // Try simple gradient decent first, then something
    // more sophisticated

    (* 
        Gradient Descent algorith:
        X_n+1 = X_n - p * ∇ f(X_n)
    *)

        /// Converts 2darray to sequence
    let arrToSeq arr =
        seq {for x in [0 .. (Array2D.length1 arr ) - 1] do 
                seq { for y in [0 .. (Array2D.length2 arr ) - 1] do
                        yield arr.[x,y] }
            }

    let multiplyMatWithConstant alpha mat  = 
        // multiplies gradient matrix by constant "alpha"
        let mapRow (row : float array) = 
            Array.map (fun x -> x * alpha) row
        in
        Array.map mapRow mat

    // let subtractMatrices mat1 mat2 = 
    //     let addRows (row1 : float array) (row2 : float array) : float array = 
    //         Array.map2 ( - ) row1 row2   
    //     in 
    //     Array2D.map2 (- ) mat1 mat2 

    let getNext (step : double) coords grad = 
        let (Gradient g) = grad in
        let (Coordinates c) = coords in
        let cmat = matrix (arrToSeq c) 
        let gradmat = matrix g
        let result = (cmat - (step * gradmat))
        Coordinates (result.ToArray ())

    let getMaxGradient (grad : Gradient) = 
        // largest absolute element of the Gradient array
        let (Gradient g) = grad
        let getMaxOfRow row = Array.max (Array.map abs row ) in 
        g |> Array.map getMaxOfRow |> Array.max


    let checkConvergence grad tol = 
        (getMaxGradient grad) < tol

    /// Turns array of coordinates into input text. 
    let createTextCoords (c : Coordinates) (txtCord : Coord) = 
        let (Coordinates cordFloats) = c in 
        let cordFloatsList = cordFloats |> arrToSeq |> Seq.map Seq.toList |> Seq.toList 
        let (Coord cordTxt) = txtCord in
        let cordTxtSplit = 
            cordTxt 
            |> List.map (fun x -> x.Split " " |> Array.filter (fun y -> y <> "") |> Array.toList) in
        let replacetxtrow (rowTxt : string list) rowFloats = 
            rowTxt.[0..1] @ (List.map string rowFloats) |> List.map (fun x -> x + "  ")
        in
        let newTxtCoord = List.map2 replacetxtrow cordTxtSplit cordFloatsList in 
        newTxtCoord |> List.map (fun x -> List.reduce (fun y z -> y + z) x) |> Coord

    let printTextCoord (c : Coordinates) (txtCoord : Coord) =
        let (Coord cList ) = createTextCoords c txtCoord in 
        List.reduce ( fun x y -> x + " \n" + y) cList

module NewtonRaphson = 
    // implement newton raphson update of geometry

    open Output.ReadOutput
    open Input.InputMBE
    open GradientDescent
    // open FSharp.Stats 
    // open Extreme.Mathematics.Optimization
    // open MathNet.Numerics.Optimization 
    open MathNet.Numerics.LinearAlgebra


    /// Compute next coordinates using NR method
    /// grad: Current step's gradient 
    /// hessInv : approx. inverse hessian
    /// coords: current step's coordinates 
    let newtonGetNext grad (hessInv : double Matrix) coords = 
        let (Gradient g) = grad in
        let (Coordinates c) = coords in
        let l1 = Array2D.length1 c in
        let l2 = Array2D.length2 c in
        let gmat = g |> Array.concat |> (fun x -> [|x|]) |> matrix
        let coordMatrix = matrix (c |> arrToSeq) 
        let deltaXflat = -1.0 * (hessInv * gmat.Transpose())
        let deltaX = Array2D.init l1 l2 (fun i j -> deltaXflat.[ i*l2 + j,0]) 
                    |> Matrix.Build.DenseOfArray
        let result = coordMatrix + deltaX
        Coordinates (result.ToArray())
