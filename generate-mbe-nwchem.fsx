let args = fsi.CommandLineArgs

let file  = args.[1]
let nbodies = int64 args.[2]
let n1bodies = int64 args.[3]

type Coord = Coord of string list
type Coordinates = Coordinates of float array array

type MBECoord = MBECoord of Coord list 

type TxtInput = TxtInput of string list 

type TotalInput = 
    {
        preamble : TxtInput
        coords :  Coord 
        postamble : TxtInput
    }
let findText (text : string) (substring : string) = 
    // returns true if "substring" exist in line of "text"
    let ind = substring.IndexOf text
    ind <> -1

let parse_coords (inputlines : string list) =
    // separates coordinates from rest of file
    let dataline = List.findIndex (findText "geometry") inputlines in
    let endline = List.findIndex (findText "end") inputlines.[dataline ..] in
    let rtrn = 
        {preamble=TxtInput inputlines.[ .. dataline];
        coords=Coord inputlines.[(dataline+1) .. (dataline + endline - 1)];
        postamble=TxtInput inputlines.[(dataline + endline) .. ]} 
    in
    rtrn

/// Change the "Start" key for nwchem inputs
let changeStart (inp : TotalInput) title = 
    let (TxtInput preamb) = inp.preamble in 
    // asume "start" is the first thing in file
    let startline = List.findIndex (findText "start") preamb in 
    let newStart = ("start " + title) ::  preamb.[(startline + 1) ..] in 
    {inp with preamble= (TxtInput newStart)}

let convertToCoord textcords = 
    // convert text coordinates to float array array
    let (Coord strlst) = textcords 
    let arrcoord = strlst |> List.toArray 
    let splitarr = arrcoord |> Array.map (fun (x : string) -> (x.Split " ") )  
                            |> Array.map (Array.filter (fun x -> x <> "") ) in 
    let c = splitarr |> Array.map (fun x -> x.[2 .. 4]) |> Array.map (Array.map float) in
    Coordinates c

let factorial (n : int64) =
    // helper function
    let rec aux (n : int64) acc = 
        match n with
        | 0L -> acc 
        | _ -> aux  (n- (int64 1)) (acc*n)
    in 
    aux n 1

let combinations n r = 
    // combinations of n objects in r groups
    (factorial n) / ( (factorial r) * (factorial (n - r)) )

let nListSelect (n : int64) (n1 : int64) (inputs : 'a list) =
    (*
        select up to n-body files
        n1 = number of 1-bodies
        inputs = list of MBE generated files
    *)
    if n>n1 then failwith "Impossible n-body selection" else
    let rec aux (n : int64) (n1 : int64) : int64 =
        match n with
        | 0L -> 0
        | x -> (combinations n1 x) + aux (n -  1L) n1
    in 
    inputs.[0 .. (int (aux n n1) - 1)]

let sortedPowerset cd = 
    // creates powerset of any list and sorts it
    let rec pwst cd = 
        match cd with
        | [] -> [[]] 
        | hd :: tl -> (List.map ((@) [hd]) (pwst tl)) @ (pwst tl)
    in pwst cd |> List.sortBy (fun s -> List.length s) |> List.tail

let limitedPowerSet (n :int64) (n1 : int64) cd=
    cd |> sortedPowerset |> nListSelect n n1

let reversePowerset cd = 
    sortedPowerset cd |> List.rev

let smallPowerset cd = 
    // creates powerset without largest term
    let rec pwst cd = 
        match cd with
        | [] -> [[]] 
        | hd :: tl -> (List.map ((@) [hd]) (pwst tl)) @ (pwst tl)
    in pwst cd |> List.tail |> List.sortBy (fun s -> List.length s) |> List.tail

let add_to_end (coords : Coord)  = 
    let (Coord x) = coords in 
    let rec aux lst n = 
        match lst with
        | [] -> []
        | x::xs -> x + "  #  " + (string n) :: aux xs (n+1)
    in Coord (aux x 1)

let generate_mbe_coords (coords : Coord) n n1  = 
    let (Coord x) = coords in 
    MBECoord ((limitedPowerSet n n1 x) |> List.map Coord)

let createFileList (mbecoords : MBECoord) (file :string) =
    let atomnumber (str : string) = 
        // returns the number after "#" in the coords
        let ind = str.IndexOf("#")
        (int str.[(ind+1) ..])
    in 
    let rec coord_num_ext coords (acc : string) = 
        // construts string from bodies that are contained in file
        match coords with 
        | [] -> acc 
        | x :: xs -> coord_num_ext xs (acc + "_" + (string (atomnumber x) + "_"))
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
    // concatenes extensions with base file name (i.e. without the extension)
        match ext with
        | [] -> []
        | x :: xs -> (f.[0 .. (f.IndexOf ".nw" - 1)] + x + "mbe.nw") :: (createFileNames xs f)
    in
    let (MBECoord x) = mbecoords
    let extensions = createExtension x in   
    createFileNames extensions file 

let createAtomArray (mbecoords : MBECoord) =
    // creates array of arrays containing atom indexes
    // eg: [[1;2];[1;3];...[1,2,3]....] 
    let atomnumber (str : string) = 
        // returns the number after "#" in the coords
        let ind = str.IndexOf("#")
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


let lines = System.IO.File.ReadLines(file) |> Seq.toList
let (pre, coor, post) = lines |> parse_coords 
                        |> (fun x -> (x.preamble, (add_to_end x.coords), x.postamble)) 

let mbe_coords = generate_mbe_coords coor nbodies n1bodies
let mbe_files = add_desc_to_coord pre post mbe_coords
let filelist = createFileList mbe_coords file

let mbe_files_corrected = List.map2 changeStart mbe_files filelist

write_to_files filelist mbe_files_corrected