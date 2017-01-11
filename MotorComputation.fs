//I'm on mono and I am limited in term of library, so let's do matrix computation manually...

//Well actually, it took more time than expected, so I decided to switch back to Python. 
//This is interesting to make in F# though, so I will probably come back to this later

///<Summary>Move an inertia matrix from the gravity center to another point</Summary>
///<param = Ig>Inertia matrix at the gravity center</param>
///<param = m>Mass of the solid</param>
///<param = a>x coordinate of the OG vector</param>
///<param = a>y coordinate of the OG vector</param>
///<param = a>z coordinate of the OG vector</param>
let HuygensMove Ig m a b c =
    let mapping i j =
        match (i,j) with
        | (0,0) -> fun x -> x + m * (b**2.0 + c**2.0)
        | (1,1) -> fun x -> x + m * (c**2.0 + a**2.0)
        | (2,2) -> fun x -> x + m * (a**2.0 + b**2.0)
        | (2,1) |(1,2) -> fun x -> x - m*b*c
        | (0,2) |(2,0) -> fun x -> x - m*c*a
        | (0,1) |(1,0) -> fun x -> x - m*a*b
    Array2D.mapi mapping Ig

let AddMatrix (M1:float[,]) (M2:float[,]) =
    Array2D.init 3 3 (fun i j -> M1.[i,j] + M2.[i,j])

let SubtractMatrix (M1:float[,]) (M2:float[,]) =
    //Reversed to facilitate pipe
    Array2D.init 3 3 (fun i j -> M2.[i,j] - M1.[i,j])

let MatrixMultiply (M1:float[,]) (M2:float[,]) =
    let getcomponent i k =
        let rec loop dim j sum =
            match j with
            | n when n = dim -> 
                sum + M1.[i,j] * M2.[j,k]
            | _ -> loop dim (j+1) (sum + M1.[i,j] * M2.[j,k])
        loop ((Array2D.length2 M1) - 1) 0 0.0
    Array2D.init (Array2D.length1 M1) (Array2D.length2 M2) getcomponent 

let MatrixTranspose (M:float[,]) = 
    Array2D.init (Array2D.length1 M) (Array2D.length2 M) (fun i j -> M.[j,i])

let MatrixDeterminant (M:float[,]) =
    let rec looprow c r nrow det =
        match r with
        | nrow -> det + M.[r,c]
        | _ -> looprow c (r+1) nrow (det + M.[r,c])
    let rec loopcol c nrow ncol det =
        match c with
        | ncol -> det + looprow c 0 nrow det
        | _ -> loopcol (c+1) nrow ncol (det + looprow c 0 nrow det)
    loopcol 0 (Array2D.length1 M) (Array2D.length2 M) 0.0

let MatrixTimesConstant (M:float[,]) (c:float) =
    Array2D.init (Array2D.length1 M) (Array2D.length2 M) (fun i j -> c * M.[i,j])

let GetCofactor (M:float[,]) = 
    let GetTempMatrix i j ic jc=
        match i with
        | n when n >= ic ->
            match j with
                |N when N >= jc -> M.[i+1,j+1]
                | _ -> M.[i+1,j]
        | _ -> 
            match j with
                | N when N >= jc -> M.[i,j+1]
                | _ -> M.[i,j]

    let cofac I J =
        let TempM = Array2D.init ((Array2D.length1 M) - 1) ((Array2D.length2 M) - 1) (fun i j -> GetTempMatrix i j I J )
        MatrixDeterminant TempM

    Array2D.init (Array2D.length1 M) (Array2D.length2 M) cofac

let MatrixInverse (M:float[,]) = 
    let det = MatrixDeterminant M
    (1.0/det) |> (MatrixTimesConstant (M |> GetCofactor |> MatrixTranspose)) 


let MakeBaseChangeMatrix Xx Xy Xz Yx Yy Yz Zx Zy Zz =
    //Cap letter = old base
    //small letter = new base
    let attributeval i j =
        match (i,j) with
        | (0,0) -> Xx
        | (0,1) -> Yx
        | (0,2) -> Zx
        | (1,0) -> Xy
        | (1,1) -> Yy
        | (1,2) -> Zy
        | (2,0) -> Xz
        | (2,1) -> Yz
        | (2,2) -> Zz
    Array2D.init 3 3 attributeval

let CylinderInertiaMatrixAtG m r h =
    let A = m * ((r**2.0)/4.0 + (h**3.0)/12.0)
    let C = m * (r ** 2.0)/2.0
    let selectValue i j =
        match i with
        |n when n = j -> 
            match n with
            | N when N = 2 -> C
            | _ -> A
        |_ -> 0.0
    Array2D.init 3 3 selectValue

let CylinderInertiaMatrixAtO m r h a b c =
    let Ig = CylinderInertiaMatrixAtG m r h
    HuygensMove Ig m a b c

let tempC r h a b c = 
    CylinderInertiaMatrixAtO (7800.0 * h/1000.0 * System.Math.PI * ((r/1000.0) ** 2.0) ) (r/1000.0) (h/1000.0) (a/1000.0) (b/1000.0) (c/1000.0)
let C1 = tempC 25.0 50.0 0.0 0.0 25.0
let C2 = tempC 60.0 15.0 0.0 0.0 -7.5
let C3 = tempC 7.0 15.0 44.0 0.0 -7.5
let C4 = tempC 7.0 15.0 (44.0 * System.Math.Cos(2.0 * System.Math.PI / 3.0)) (44.0 * System.Math.Sin(2.0 * System.Math.PI / 3.0)) -7.5
let C5 = tempC 7.0 15.0 (44.0 * System.Math.Cos(4.0 * System.Math.PI / 3.0)) (44.0 * System.Math.Sin(4.0 * System.Math.PI / 3.0)) -7.5
let C6 = tempC 7.0 30.0 44.0 0.0 -15.0

let C = C1 |> (AddMatrix C2) |> (SubtractMatrix C3) |> (SubtractMatrix C4) |> (SubtractMatrix C5) |> (AddMatrix C6)


let printmatrix (M:float[,]) =
    let Mf = Array2D.init (Array2D.length1 M) (Array2D.length2 M) (fun i j -> sprintf "%f" (1000000.0 * M.[i,j]))
    printf "%A" Mf
    printf "\n"

//printf "%A" (C)
printmatrix C

let BaseChange = MakeBaseChangeMatrix 0.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0
printmatrix BaseChange
printmatrix (MatrixMultiply C BaseChange)
((Array2D.init 3 1 (fun i j -> 
    match i with
    | 2 -> 1.0
    | _ -> 0.0)) |> (MatrixMultiply BaseChange)) |> printmatrix

printmatrix (GetCofactor C)

printmatrix ((MatrixInverse C) |> (MatrixMultiply C))