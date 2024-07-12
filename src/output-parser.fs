namespace Output

module ReadOutput =
    open External.ExecuteExternal
    open Input.InputMBE
    open MathNet.Numerics.LinearAlgebra 

    type GmsOutput = GmsOutput of string array

    type Gradient = Gradient of float array array

    type Energy = Energy of float


    let gmsOutput (outputRec : CommandResult) = 
        GmsOutput (outputRec.StandardOutput.Split "\n")

    let printOutput output file =
        Array.iter (fprintfn file "%s") output

    let getGrad (output : GmsOutput) = 
        let errorFile = "error.log"
        let (GmsOutput out) = output  in 
        let gradHead = 
            try Array.findIndex (findText "UNITS ARE HARTREE/BOHR") out
            with 
            | :? System.Collections.Generic.KeyNotFoundException -> begin
                    use er = System.IO.File.AppendText errorFile
                    fprintfn er "Exception, phrase <UNITS ARE HARTREE/BOHR> not found"
                    fprintfn er " File is \n"
                    printOutput out er
                    0 
                end
        in
        let gradEnd = 
            try Array.findIndex (findText "MAXIMUM GRADIENT =") out 
            with 
            | :? System.Collections.Generic.KeyNotFoundException -> begin
                    use er = System.IO.File.AppendText errorFile
                    fprintfn er "Exception, phrase <MAXIMUM GRADIENT> not found"
                    fprintfn er " File is \n"
                    printOutput out er  
                    0 
                end
        in
        let gradRaw = out.[(gradHead + 1) .. (gradEnd - 2)] in 
        let splitarr = gradRaw |> Array.map (fun (x : string) -> (x.Split " ") )  
                               |> Array.map (Array.filter (fun x -> x <> "") ) in 
        let gradient = splitarr |> Array.map (fun x -> x.[2 ..]) |> Array.map (Array.map float) |> Gradient in
        gradient

    let getEnergy (output : GmsOutput) = 
        // For now only searches MP2 energy. 
        let (GmsOutput out) = output in 
        // go over twice to make sure you get the right "Total Energy" string
        let head1 = Array.findIndex (findText "MP2 PROPERTIES") out in
        let head2 = Array.findIndex (findText "TOTAL ENERGY") out.[head1 ..] in
        let energyLine = out.[head1 + head2] in 
        let splitline = energyLine.Split " " |> Array.filter (fun x -> x <> "") in
        let energy = splitline.[splitline.Length - 1] |> float in energy 

    /// Make a gradient of an n-body term equal in size to total gradient. 
    let equalGradDims (gradient : Gradient) (index : int list) (n : int) = 
        let (Gradient grad) = gradient in 
        let initgrad = [| for i in 1 .. n do Array.create 3 0.0 |] in
        let mutable gradi = 0 in
        for i in index do 
            initgrad.[i-1] <- grad.[gradi]
            gradi <- gradi + 1
        Gradient initgrad 


    let addGradientsMatrices (mat1 : Gradient) (mat2 : Gradient) =
    // add 2 matrices. Must be 2-dimensional
        let (Gradient grad1) = mat1 in
        let (Gradient grad2) = mat2 in
        let addRows (row1 : float array) (row2 : float array) : float array = 
            Array.map2 ( + ) row1 row2   
        in 
        Array.map2 addRows grad1 grad2 |> Gradient  

    let subtractGradientsMatrices (mat1 : Gradient) (mat2 : Gradient) =
    // add 2 matrices. Must be 2-dimensional
        let (Gradient grad1) = mat1 in
        let (Gradient grad2) = mat2 in
        let subtractRows (row1 : float array) (row2 : float array) : float array = 
            Array.map2 ( - ) row1 row2   
        in 
        Array.map2 subtractRows grad1 grad2 |> Gradient  
    
    let gradientDifferenceToMatrix (grad1 : Gradient) (grad2 : Gradient) = 
            let (Gradient temp) = subtractGradientsMatrices grad1 grad2 in 
            let x = matrix temp
            Matrix.map double x

