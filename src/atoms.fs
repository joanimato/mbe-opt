namespace Atoms
(* 
     COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
     AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT
     HE,NE,AR,KR,RN ARE FROM: "COVALENT RADII REVISITED".
     BEATRIZ CORDERO, VERONICA GOMEZ, ANA E. PLATERO-PRATS,
     MARC REVES, JORGE ECHEVERRIA, EDUARD CREMADES,
     FLAVIA BARRAGAN, SANTIAGO ALVAREZ
     DALTON TRANS.: 2832-2838(2008)   DOI:10.1039/B801115J 
     
     and

     "MOLECULAR SINGLE-BOND COVALENT RADII FOR ELEMENTS 1-118".
     P. PYYKKO, M. ATSUMI
     CHEMISTRY: A EUROPEAN JOURNAL 15: 186-197(2009)
     DOI:10.1002/CHEM.200800987.

     *)

module AtomDefs = 
    open Input.InputMBE

    type Atom = {
        name: string
        symbol: string
        mass : double  // Atomic mass
        znum : int  // Atomic number
        vdwr :  double option // van der Waals radii
        rcov : double option // Covalent radii
    }

    type AtomCoordinates = {
        atom : Atom
        xcoord : double
        ycoord : double
        zcoord : double
    } 

    type Molecule = Molecule of AtomCoordinates array
    // wdwr in Angstroms

    let hydrogen = 
        {name = "Hydrogen"; symbol = "H"; mass=1.008; znum=1; vdwr=Some 1.2;
        rcov= Some 0.3}
    let helium =
        {name = "Helium"; symbol = "He"; mass=4.002602; znum=2; vdwr=Some 1.4 ;
        rcov=Some 0.46}
    let lithium =
        {name = "Lithium"; symbol = "Li"; mass=6.94; znum=3; vdwr=Some 1.82 ;
        rcov=Some 1.23}
    let beryllium =
        {name = "Beryllium"; symbol = "Be"; mass=9.012183; znum=4; vdwr=Some 1.53 ;
        rcov=Some 0.89}
    let boron =
        {name = "Boron"; symbol = "B"; mass=10.81; znum=5; vdwr=Some 1.92 ;
        rcov=Some 0.88}
    let carbon =
        {name = "Carbon"; symbol = "C"; mass=12.011; znum=6; vdwr=Some 1.7 ;
        rcov=Some 0.77}
    let nitrogen =
        {name = "Nitrogen"; symbol = "N"; mass=14.007; znum=7; vdwr=Some 1.55 ;
        rcov= Some 0.7}
    let oxygen =
        {name = "Oxygen"; symbol = "O"; mass=15.999; znum=8; vdwr=Some 1.52 ;
        rcov=Some 0.66}
    let fluorine =
        {name = "Fluorine"; symbol = "F"; mass=18.998; znum=9; vdwr=Some 1.47 ;
        rcov=Some 0.58}
    let neon =
        {name = "Neon"; symbol = "Ne"; mass=20.1797; znum=10; vdwr=Some 1.54 ;
        rcov=Some 0.67}
    let sodium =
        {name = "Sodium"; symbol = "Na"; mass=22.989; znum=11; vdwr=Some 2.27 ;
        rcov=Some 1.66}
    let magnesium =
        {name = "Magnesium"; symbol = "Mg"; mass=24.305; znum=12; vdwr=Some 1.73 ;   
        rcov=Some 1.36 } 
    let potassium =
        {name = "Potassium"; symbol = "K"; mass=39.0983; znum=19; vdwr=Some 2.75 ;
        rcov=Some 2.03}
    let calcium =
        {name = "Calcium"; symbol = "Ca"; mass=40.078; znum=20; vdwr=Some 2.31 ;
        rcov=Some 1.74}
    let zinc =
        {name = "Zinc"; symbol = "Zn"; mass = 65.38; znum=30; vdwr = Some 1.39;
        rcov=Some 1.25}
// FOR BADGER'S RULES SEE J.C.P. 2, 128(1934)

// A GENERALIZED BADGER'S RULE:  D.R.HERSCHBACH, V.W.LAURIE, J.CHEM.PHYS. 35, 458-463(1961).

    ///match individual arrays of coordinates with predifned atoms. 
    let matchAtoms (atArray : string array) = 
        let createAtomCoordinate atomname = 
            let doubleCoords = Array.map double atArray.[2..4] in
            {atom=atomname; xcoord=doubleCoords.[0];
            ycoord=doubleCoords.[1];zcoord=doubleCoords.[2]}
        match atArray.[0].ToLower() with 
        | "h" -> createAtomCoordinate hydrogen 
        | "he" -> createAtomCoordinate helium
        | "li" -> createAtomCoordinate lithium
        | "be" -> createAtomCoordinate beryllium
        | "na" -> createAtomCoordinate sodium
        | "mg" -> createAtomCoordinate magnesium
        | "k" -> createAtomCoordinate potassium
        | "ca" -> createAtomCoordinate calcium
        | "zn" -> createAtomCoordinate zinc
        | x -> failwith $"unknown atom type %A{x}"


    /// Convert Coordinate Array into a Molecule type
    let coordToMolecule (coord : Coord) = 
        let (Coord c) = coord
        let arrcoord = c |> List.toArray 
        let splitarr = arrcoord |> Array.map (fun (x : string) -> (x.Split " ") )  
                               |> Array.map (Array.filter (fun x -> x <> "") ) in 
        let molecule = Array.map matchAtoms splitarr
        Molecule molecule

module Hessian = 

    open External.ExecuteExternal
    open AtomDefs
    open System 
    open MathNet.Numerics.LinearAlgebra 

    // do a semiempirical calculation to get external hessian

    /// Guess hessian from empirical parameters 
    let hessGuessCartesian (mol : Molecule) = 
        let (Molecule molecule) = mol in
        let nAtoms = Array.length molecule in 
        // populate initial 3N*3N hessian
        let hess = 
            Array2D.create (nAtoms*3) (nAtoms*3) (double 0.0) in

        // tempporary storage 
        let mutable temp =  Array2D.create 3 3 (double 0.0)

        // loop variables
        let mutable rx = double 0.0
        let mutable ry = double 0.0 
        let mutable rz = double 0.0
        let mutable rr = double 0.0
        let mutable rab = double 0.0 
        let mutable cab = double 0.0 

        // some indices 
        let mutable i2 = 0
        let mutable j2 = 0

        for i = 0 to nAtoms - 1 do 
            for j = 0 to (i - 1) do


                rx <- molecule.[i].xcoord - molecule.[j].xcoord
                ry <- molecule.[i].ycoord - molecule.[j].ycoord
                rz <- molecule.[i].zcoord - molecule.[j].zcoord
                rr <- rx**2.0 + ry**2.0 + rz**2.0
                rab <- Math.Sqrt rr
                cab <- molecule.[i].atom.rcov.Value + molecule.[i].atom.rcov.Value
                // hessian with respect to rab
                let hess2 = (double 0.3601) * Math.Exp ((double -1.944) * (rab-cab) ) in
                // convert to cartesian with "approximate" chain rule
                // from gamess , I think corners are cut here
                temp.[0,0] <- (rx*rx*hess2)/rr
                temp.[0,1] <- (rx*ry*hess2)/rr
                temp.[0,2] <- (rx*rz*hess2)/rr
                temp.[1,1] <- (ry*ry*hess2)/rr
                temp.[1,2] <- (ry*rz*hess2)/rr
                temp.[2,2] <- (rz*rz*hess2)/rr
                temp.[1,0] <- temp.[0,1]
                temp.[2,0] <- temp.[0,2]
                temp.[2,1] <- temp.[1,2]

                i2 <- (i )*3
                j2 <- (j )*3
                //printfn "i2: %d ; j2: %d" i2 j2
                // construct xyz hessian
                for ii = 0 to 2 do
                    for jj = 0 to 2 do 
                        hess.[i2+ii,i2+jj] <- (hess.[i2+ii,i2+jj] + temp.[ii,jj])
                        hess.[j2+ii,j2+jj] <- (hess.[j2+ii,j2+jj] + temp.[ii,jj])
                        hess.[i2+ii,j2+jj] <- (hess.[i2+ii,j2+jj] - temp.[ii,jj])
                        hess.[j2+ii,i2+jj] <- (hess.[j2+ii,i2+jj] - temp.[ii,jj])
                 
        // check for small diagonal elements
        for i=0 to (nAtoms*3 - 1) do
            if (hess.[i,i] < 1.0e-5 ) then hess.[i,i] <- double 1.0e-5

        Matrix.Build.DenseOfArray hess


    /// Updates the inverse of Hessian using BFGS algorithm.
    /// hessInv -> Inverted hessian
    /// deltaGrad -> Difference between new and old gradient
    /// deltaX -> difference between new and old coordinates
    let hessianInverseUpdateBFGS (hessInv : Matrix<double>) (deltaG : double Matrix) (deltaX : double Matrix)= 
        let deltaXflat = deltaX.ToRowArrays() |> Array.concat |> (fun x -> [|x|]) |> matrix
        let deltaXTranspose= deltaXflat.Transpose() in 
        let deltaGflat = deltaG.ToRowArrays() |> Array.concat |> (fun x -> [|x|]) |> matrix
        let deltaGTranspose = deltaGflat.Transpose() in 
        let firstDenom = (deltaGflat * deltaXTranspose).[0,0] in 
        let firstTerm = deltaXTranspose * deltaXflat / firstDenom in 
        let secondDenom : double = (deltaGflat * hessInv * deltaGTranspose).[0,0] in 
        let secondTerm = 
            (hessInv * deltaGTranspose * deltaGflat * hessInv) / secondDenom in 
        let newHessian = hessInv + (firstTerm  - secondTerm)
        newHessian

    
        
