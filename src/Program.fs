open Input.InputMBE
open Output.ReadOutput
open External.ExecuteExternal
open Optimize.GradientDescent
open Optimize.NewtonRaphson
open Argu
open Atoms
open MathNet.Numerics.LinearAlgebra
// open MathNet.Numerics.Probability


// Create CLI argument parser
type CliArguments = 
    | Input_File of string
    | Nbody of int 
    | Gms_Exec of string 
    | Iter of int 
    | Ncpus of int

    interface IArgParserTemplate with
        member s.Usage = 
            match s with 
                | Input_File _ -> "specify input file"
                | Nbody _ -> "highest n-body to include"
                | Gms_Exec _ -> "location of gamess executable"
                | Iter _ -> "max iterations"
                | Ncpus _ -> "cpus per calculation"


type Inputfile = Inputfile of string 

type InputMbeFiles = InputMbeFiles of Inputfile list

let count1Bodies mbecoords =
    // counts how many 1-bodies exist. Assumes sorted list. 
    let (MBECoord x) = mbecoords in 
    let rec aux acc mcords = 
        match mcords with
        | [] -> acc
        | x :: xs ->
            begin
                let (Coord xlist) = x
                if (List.length xlist) = 1 then (aux (acc+1) xs) else acc
            end
    in
    aux 0 x


let nBodySelect n n1 (inputs : InputMbeFiles) =
    (*
        select up to n-body files
        n1 = number of 1-bodies
        inputs = list of MBE generated files
    *)
    if n>n1 then failwith "Impossible n-body selection" else
    let rec aux n n1 =
        match n with
        | 0 -> 0
        | x -> (combinations n1 x) + aux (n-1) n1
    in 
    let (InputMbeFiles inputlist) = inputs
    InputMbeFiles inputlist.[0 .. ((aux n n1) - 1)]

let executeGamessAsync (gms : string) (files : InputMbeFiles) ncpus = 
    let (InputMbeFiles nbodies) = files
    let n = string ncpus
    let runs = seq {for f in nbodies do 
                        let (Inputfile fstring) = f;
                        executeCommand gms [fstring; n] } 
    let results = runs |> Async.Parallel |> Async.RunSynchronously in
    results


let combineGradients grads_with_index indeces n =
    // calculates gradient of 1 term in the mbe
    // Requires:
    // grads_with_index = List of tules containing gradients and index of term
    // indeces = just the indeces of nbody terms
    // n = size of the body. 
    let zeroMatrix = Gradient (Array.create n (Array.create 3 0.0))
    let findgradient (ilst : int list) (gwi : (int list * Gradient) list )  = 
        let (_, grad) = List.find (fun (x,_) -> x = ilst) gwi in
        grad 
    in
    let rec aux2 ind1 index gwi =
        // assume largest element first 
        // ind1 is the first element of index
        // gwi stays the same. 
        match index with 
        | [] -> zeroMatrix
        | x :: xs when (List.length x) = 1 -> zeroMatrix
        | x :: xs -> 
            begin
                let powerx = smallPowerset x in
                subtractGradientsMatrices (findgradient ind1 gwi) 
                    (List.fold addGradientsMatrices zeroMatrix 
                        (List.map (fun z -> aux2 z (List.rev (sortedPowerset z)) gwi) powerx)) 
            end
    in 
    let ind1 = List.head indeces in 
    aux2 ind1 indeces grads_with_index 

let combineEnergies ens_with_indices indeces =
    let zeroEn = 0.0 in
    let findEnergy (ilst : int list) (ewi : (int list * float) list )  = 
        let (_, en) = List.find (fun (x,_) -> x = ilst) ewi in
        en
    in
    let rec aux ind1 index ewi  =
        // assume largest element first 
        // ind1 is the first element of index
        // gwi stays the same. 
        match index with 
        | [] -> zeroEn
        | x :: xs when (List.length x) = 1 -> findEnergy x ewi
        | x :: xs ->
            begin
                let powerx = smallPowerset x in
                    (findEnergy ind1 ewi) - 
                    (List.fold ( + ) zeroEn
                    (List.map (fun z -> aux z (List.rev (sortedPowerset z)) ewi) powerx)) 
            end
    in 
    aux (List.head indeces ) indeces ens_with_indices
  


let createLogFile logfile =
    use outFile = System.IO.File.CreateText logfile
    fprintfn outFile "Starting MBE calculation"

let updateLogFile logfile txt = 
    use outFile = System.IO.File.AppendText logfile
    fprintfn outFile txt

let combineEnergiesFormula s = failwith "Implement"

let combineGradientsFormula s = failwith "Implement"

/// Obtain MBE energy out of formatted output
/// mbeIndeces -> array array of mbe indices
/// nbody -> nbody to compute the energy
/// n1body -> number of 1-bodies
let getTotalEnergy mbeIndeces formatedout nbody n1body = 
    let indexlist = mbeIndeces |> Array.toList |> List.map (fun x -> Array.toList x) |> nListSelect nbody n1body
    let revIndexList =  indexlist |> List.rev 
    let energies  = Array.map getEnergy formatedout |> Array.toList |> List.rev
    let ens_with_indices  = List.map2 (fun x y -> (x,y) ) revIndexList energies
 
    let mutable totalEnergy = 0.0 

    for jl in revIndexList do
        totalEnergy <- (combineEnergies ens_with_indices (reversePowerset jl)) + totalEnergy
    totalEnergy 

/// take next step NR in optimization 
let getTotalGrad gradients mbeIndeces n1body nbody =
    let indexlist = mbeIndeces |> Array.toList |> List.map (fun x -> Array.toList x) |> nListSelect nbody n1body
    let revIndexList =  indexlist |> List.rev 
    let sameDimGradients = 
            List.map2 (fun x  y -> equalGradDims x y n1body) gradients indexlist 
            |> List.rev
    let revGradients = gradients |> List.rev
    let grads_with_indeces = List.map2 (fun x y -> (x,y)) revIndexList sameDimGradients  
  
    (* ************************
        Combine gradients
    ************************** *)  
    let mutable totalgrad = Gradient (Array.create n1body (Array.create 3 0.0)) in

    for jl in revIndexList do
        totalgrad <- (combineGradients grads_with_indeces (reversePowerset jl) n1body ) 
                        |> addGradientsMatrices totalgrad

    //let newCoords = newtonGetNext totalgrad hess coords

    totalgrad

[<EntryPoint>]
let main argv = 

    // Parse command line arguments with Argu
    let parser = ArgumentParser.Create<CliArguments>(programName = "mbe-opt")
    let usage = parser.PrintUsage

    //failwith "got here"
    let results = parser.ParseCommandLine (inputs=argv, raiseOnUsage = true) 
    let file = results.GetResult Input_File
    let nbody = results.GetResult Nbody 
    let gamessExec = results.GetResult Gms_Exec
    let itermax = results.GetResult Iter 
    let ncpus = results.GetResult Ncpus
    //let allResults = results.GetAllResults

    // let file = argv.[0]  // initial input file
    // let nbody = int argv.[1] // n-body term to do expansion calculation
    // let gamessExec = argv.[2] // complete path of gamess executable 
    // let itermax = int argv.[3] // maximum number of geom. iterations

    // temporary testing
    // let file = "li4-mul1-mp2-lin.inp"
    // let file = "exam04.inp"
    // let file = "be4-mp2-td.inp"
    // let nbody = 2
    // let gamessExec = "/Users/mato166/gamess-dev/gms.sh"
    // let itermax = 100 
    // let file = "be4-mul1-dft-lin.inp"

    (* ************************************
        Create MBE of initial coordinates 
    ************************************* *)

    // Write output in output.log
    
    let outFile = "output.log"
    let errorFile = "error.log"
    createLogFile outFile 
    createLogFile errorFile

    let lines = System.IO.File.ReadLines(file) |> Seq.toList
    let (pre, coor, post) = lines |> parse_coords 
                            |> (fun x -> (x.preamble, (add_to_end x.coords), x.postamble)) 

    let mol = AtomDefs.coordToMolecule coor

    updateLogFile outFile $"File read: %s{file} \n\n"
    updateLogFile outFile "Coordinates extracted from main file \n"

    let n1body =  // temporary --> treat all atoms as 1bodies
        let (AtomDefs.Molecule c) = mol in Array.length c

    let mbe_coords = generate_mbe_coords coor nbody n1body
    let mbe_files = add_desc_to_coord pre post mbe_coords
    let filelist = createFileList mbe_coords file
    let mbeIndeces = createAtomArray mbe_coords |> List.toArray

    write_to_files filelist mbe_files

    updateLogFile outFile  " Initial MBE files generated. \n "

    let mbeFileList = InputMbeFiles (List.map Inputfile filelist)
 

    let nbodyMBEFiles = nBodySelect nbody n1body mbeFileList

    updateLogFile outFile $"The MBE for %d{nbody} out of %d{n1body} bodies will be computed \n"

    (* ******************
      Take First Step
    ******************** *)
  
    //  some variables to keep track of //
    let mutable coordinates = convertToCoord coor 
    let mutable prevCoordinates = convertToCoord coor
    let mutable diffCoord = coordToMatrix coordinates
    let mutable prevEnergy = 0.0
    let mutable totalEnergy = 0.0 
    let mutable enDiff = 0.0 
    let mutable mainLoopBreak = true
    let mutable iter = 0
    let mutable newCoor = coor
    let mutable totalgrad = Gradient (Array.create n1body (Array.create 3 0.0)) in
    let mutable prevGrad = Gradient (Array.create n1body (Array.create 3 0.0))
    let mutable gradDiff = matrix (Array.create n1body (Array.create 3 (double 0.0)))

    (* *****************************************
        Run GAMESS calculations in parallel. 
    ****************************************** *)
    let results = executeGamessAsync gamessExec nbodyMBEFiles ncpus

    // Generate initial hessian
    let mutable hess = Hessian.hessGuessCartesian mol
    // Format output
    let formatedout = Array.map gmsOutput results 
    let gradients = Array.map getGrad formatedout |> Array.toList  // needs

    prevGrad <- getTotalGrad gradients mbeIndeces n1body nbody

    // Initial STEP if not already converged
    if (checkConvergence prevGrad 1.0e-4) then
        begin 
            updateLogFile outFile "Geometry is already converged !! \n"
            updateLogFile outFile $"%s{(printTextCoord coordinates coor ) } \n" 
            updateLogFile outFile "Converged Gradient: \n"
            updateLogFile outFile $"%A{prevGrad}"  
            mainLoopBreak <- false // never start loop
        end
    else
        begin
            let coords = coordinates
            coordinates <- newtonGetNext prevGrad hess coords
            updateLogFile outFile "NR step taken for new coordinates.\n"
        end


    prevEnergy <- getTotalEnergy mbeIndeces formatedout nbody n1body

    let gradMax = getMaxGradient prevGrad

    updateLogFile outFile $"Iteration: 0, Energy: {prevEnergy}, E. Diff: {enDiff} Max Grad: {gradMax}" 
    updateLogFile outFile "Current Geometry: \n"
    updateLogFile outFile $"%s{(printTextCoord coordinates coor )} \n"

    // Loop over iterations for geometry optimization. 
    while mainLoopBreak do
        iter <- iter + 1

        (* ********************************
            Generate new file from coords
        ********************************* *) 

        newCoor <- createTextCoords coordinates newCoor

        let mbe_coords = generate_mbe_coords newCoor nbody n1body
        let mbe_files = add_desc_to_coord pre post mbe_coords // pre and post don't change

        write_to_files filelist mbe_files

        (* ***************************************
            Run GAMESS calculations in parallel. 
        ****************************************** *)

        let results = executeGamessAsync gamessExec nbodyMBEFiles ncpus
        // results contain array of output from all gamess calculations 

        (* **************************************
            Parse output and save gradients 
        *************************************** *)        
        
        let formatedout = Array.map gmsOutput results 
        let gradients = Array.map getGrad formatedout |> Array.toList  // needs to be equal or not?

        (* ***************************
            Combine energy/gradients
        ****************************** *)  

        totalgrad <- getTotalGrad gradients mbeIndeces n1body nbody
        totalEnergy <- getTotalEnergy mbeIndeces formatedout nbody n1body



        // calculate differences between current and previous stpe
        // and update variables
        enDiff <- totalEnergy - prevEnergy
        prevEnergy <- totalEnergy
        gradDiff <- gradientDifferenceToMatrix totalgrad prevGrad
        prevGrad <- totalgrad
        diffCoord <- (coordToMatrix coordinates) - (coordToMatrix prevCoordinates)
        prevCoordinates <- coordinates

        // Hessian update step 
        hess <- Hessian.hessianInverseUpdateBFGS hess gradDiff diffCoord

        (* ******************** 
            NR optimization
        *********************** *) 

        if (checkConvergence totalgrad 1.0e-4) then
            begin 
                updateLogFile outFile "Geometry has converged !! \n"
                updateLogFile outFile $"%s{(printTextCoord coordinates coor ) } \n" 
                updateLogFile outFile "Converged Gradient: \n"
                updateLogFile outFile $"%A{totalgrad}"  
                mainLoopBreak <- false
            end
        else
            begin
                let coords = coordinates
                coordinates <- newtonGetNext totalgrad hess coordinates
            end


        let gradMax = getMaxGradient totalgrad  
        updateLogFile outFile $"Iteration: {iter}, Energy: {totalEnergy}, E. Diff: {enDiff} Max Grad: {gradMax}" 
        updateLogFile outFile "Current Geometry: \n"
        updateLogFile outFile $"%s{(printTextCoord coordinates coor )} \n"

        // step after iteration limit 
        if iter >= itermax then begin 
            updateLogFile outFile "Max Itermations Reached."
            updateLogFile outFile "Coordinates thus far: "
            updateLogFile outFile $"%s{printTextCoord coordinates coor} \n" 
            updateLogFile outFile " "
            updateLogFile outFile "Gradient thus far \n"
            updateLogFile outFile $"%A{totalgrad}" 
            mainLoopBreak <- false
            end

    0
