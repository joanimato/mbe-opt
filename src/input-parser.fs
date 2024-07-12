namespace Input

/// Read input files for the MBE run
module InputParser =    
    
    type InputText = 
        | MbeInput of string array
        | GmsInput of string array

    
/// Process GAMESS input file 
module InputMBE = 
    open MathNet.Numerics.LinearAlgebra

    type Coord = Coord of string list
    type Coordinates = Coordinates of float [,]
    type MBECoord = MBECoord of Coord list 



    type TxtInput = TxtInput of string list 

    type TotalInput = 
       {
           preamble : TxtInput
           coords :  Coord 
           postamble : TxtInput
       }
    let findText (substring : string) (text : string) = 
        // returns true if "substring" exist in line of "text"
        // not case sensitive
        let ind = text.IndexOf(substring,System.StringComparison.CurrentCultureIgnoreCase)
        ind <> -1

    let parse_coords (inputlines : string list) =
        // separates coordinates from rest of file
        let dataline = List.findIndex (findText " $DATA") inputlines in
        let endline = List.findIndex (findText " $END") inputlines.[dataline ..] in
        let rtrn = 
            {preamble=TxtInput inputlines.[ .. (dataline + 2)];
            coords=Coord inputlines.[(dataline+3) .. (dataline + endline - 1)];
            postamble=TxtInput inputlines.[(dataline + endline) .. ]} 
        in
        rtrn


    ///get atom list from coordinate. 
    let getAtoms coor = 
        let (Coord c) = coor in
        // first element of each item in coord list
        let getfirst (item : string) = item.Split " " 
                                    |> Array.filter (fun x -> x <> "") 
                                    |> (fun x y -> Array.get y x) 0 in 
        List.map getfirst c

    let convertToCoord textcords = 
        // convert text coordinates to float array array
        let (Coord strlst) = textcords 
        let arrcoord = strlst |> List.toArray 
        let splitarr = arrcoord |> Array.map (fun (x : string) -> (x.Split " ") )  
                               |> Array.map (Array.filter (fun x -> x <> "") ) in 
        let c = splitarr |> Array.map (fun x -> x.[2 .. 4]) |> Array.map (Array.map float) in
        let d = array2D c
        Coordinates d

    let parse_desc inp_file = 
        // separates description from coordinates of file
        let dataline = Array.findIndex (findText " $DATA") inp_file in 
        inp_file.[0 .. (dataline - 1)]

    let sortedPowerset cd = 
        // creates powerset of any list and sorts it
        let rec pwst cd = 
            match cd with
            | [] -> [[]] 
            | hd :: tl -> (List.map ((@) [hd]) (pwst tl)) @ (pwst tl)
        in pwst cd |> List.sortBy (fun s -> List.length s) |> List.tail

    let reversePowerset cd = 
        sortedPowerset cd |> List.rev

    let smallPowerset cd = 
        // creates powerset without largest term
        let rec pwst cd = 
            match cd with
            | [] -> [[]] 
            | hd :: tl -> (List.map ((@) [hd]) (pwst tl)) @ (pwst tl)
        in pwst cd |> List.tail |> List.sortBy (fun s -> List.length s) |> List.tail

    let factorial n =
        // helper function
        let rec aux n acc = 
            match n with
            |0 -> acc 
            | x -> aux (n-1) (acc*n)
        in 
        aux n 1

    let combinations n r = 
        // combinations of n objects in r groups
        (factorial n) / ( (factorial r) * (factorial (n - r)) )

    let nListSelect n n1 (inputs : 'a list) =
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
        inputs.[0 .. ((aux n n1) - 1)]

    let limitedPowerSet n n1 cd=
        cd |> sortedPowerset |> nListSelect n n1

    let add_to_end (coords : Coord)  = 
        let (Coord x) = coords in 
        let rec aux lst n = 
            match lst with
            | [] -> []
            | x::xs -> x + "  !  " + (string n) :: aux xs (n+1)
        in Coord (aux x 1)

    let generate_mbe_coords (coords : Coord) n n1 = 
        let (Coord x) = coords in 
        MBECoord ((limitedPowerSet n n1 x) |> List.map Coord)

    let createFileList (mbecoords : MBECoord) (file :string) =
        let atomnumber (str : string) = 
            // returns the number after "!" in the coords
            let ind = str.IndexOf("!")
            (int str.[(ind+1) ..])
        in 
        let rec coord_num_ext coords (acc : string) = 
            // construts string from bodies that are contained in file
            match coords with 
            | [] -> acc 
            | x :: xs -> coord_num_ext xs (acc + (string (atomnumber x)))
        in
        let rec createExtension (coords : Coord list) = 
            match coords with 
            // List of List of string. Ex: [[string,string,string],[string,string],...]
            // each string item is a line of coords. . 
            | [] -> []
            | x :: xs -> begin 
                let (Coord str) = x in
                (coord_num_ext str "") :: (createExtension xs)
                end
        in
        let rec createFileNames ext (f : string) =
        // concatenes extensions with base file name (i.e. without the .inp)
            match ext with
            | [] -> []
            | x :: xs -> (f.[0 .. (f.IndexOf ".inp" - 1)] + x + "mbe.inp") :: (createFileNames xs f)
        in
        let (MBECoord x) = mbecoords
        let extensions = createExtension x in   
        createFileNames extensions file 
    
    let createAtomArray (mbecoords : MBECoord) =
        // creates array of arrays containing atom indexes
        // eg: [[1;2];[1;3];...[1,2,3]....] 
        let atomnumber (str : string) = 
            // returns the number after "!" in the coords
            let ind = str.IndexOf("!")
            (int str.[(ind+1) ..])
        in
        let rec constructOneArray (coords : Coord) (acc : int list) =
            let (Coord c) = coords in
            match c with
            | [] -> acc 
            | x :: xs -> constructOneArray (Coord xs) (atomnumber x :: acc )
        in 
        let rec constructAllArrays (allcoords : Coord list) =
            match allcoords with 
            | [] -> []
            | x :: xs -> List.rev (constructOneArray x []) :: constructAllArrays xs 
        in
        let (MBECoord c) = mbecoords
        constructAllArrays c |> List.map List.toArray 



    let add_desc_to_coord (pre : TxtInput) (post: TxtInput) (coords : MBECoord) = 
        // adds file description to each mbe coord
        let (MBECoord x) = coords in
        let rec aux pre' post' coords' =
            match coords' with 
            | [] -> []
            | x :: xs -> {preamble=pre ; coords=x ; postamble=post} :: aux pre post xs 
        in
        aux pre post x

    let combineInputFile (input : TotalInput) = 
        // combine fields of TotalInput into a single list of strings
        let (TxtInput x) = input.preamble
        let (Coord y) = input.coords
        let (TxtInput z) = input.postamble
        x @ y @ z

    let print_list filename (input : TotalInput) =
        // print elements of string list into filename
        use file = System.IO.File.CreateText(filename)
        input |> combineInputFile |> List.iter (fun elem -> fprintf file "%s \n" elem)

    let write_to_files filelist mbecoords = 
        (filelist, mbecoords) ||> List.iter2 print_list 

    let coordToMatrix coord : double Matrix= 
        let (Coordinates c) = coord in 
        Matrix.Build.DenseOfArray c

    let createMBE (inpfile : string) n n1 =
        (* 
            Take a gamess input file and output all 
            possible many-body input files 
        *)
        let file = inpfile
        // temporary testing
        let file = "li4-mul1-mp2-lin.inp"
        //let file = "be4-mul1-dft-lin.inp"

        let lines = System.IO.File.ReadLines(file) |> Seq.toList
        
        let (pre, coor, post) = lines |> parse_coords 
                                |> (fun x -> (x.preamble, (add_to_end x.coords), x.postamble)) 
        let mbe_coords = generate_mbe_coords coor n n1
        let mbe_files = add_desc_to_coord pre post mbe_coords
        let filelist = createFileList mbe_coords file

        write_to_files filelist mbe_files
