SetDirectory[NotebookDirectory[]];
info := {
  Print[StyleForm[
    "=====================================================", 
    "Section", FontSize -> 14, Black]];
  Print[StyleForm["PACKAGE:", "Section", FontSize -> 14], 
   StyleForm["LINBLAD OPERATOR", "Section", FontSize -> 14, Black] ];
  Print[StyleForm["BY: Jorge Chávez-Carlos, 2024", "Section", 
    FontSize -> 12, Black]];
  Print[StyleForm[
    "=====================================================", 
    "Section", FontSize -> 14, Black]];
  Print[StyleForm["LINBLAD_1.0 package:", "Section", FontSize -> 12, 
    Black]];}; Clear[HH, dim, Cnames, lnlindblad, CCs, gs, L, Ld, \
\[Rho], H, h, LIND, LINDBLAD];
HH = Import["H.dat"][[1 ;; -1, 1]];
dim = Sqrt[Dimensions[Flatten[HH]]][[1]];
HH = Partition[HH, dim];
Cnames = FileNames["C*.dat"];
lnlindblad = Length[Cnames];
CCs = Table[
   Partition[Import[Cnames[[i]]][[1 ;; -1, 1]], dim], {i, 1, 
    lnlindblad}];
gs = Import["g.dat"][[1 ;; -1, 1]];
Conm[A_, B_] := A . B - B . A;
AntConm[A_, B_] := A . B + B . A;
Do[
  Do[
   Do[
     L[i, j, k] = Subscript[L, i, j, k];
     Ld[i, j, k] = Subscript[L, i, j, k];
     , {j, 0, dim - 1, 1}];
   , {i, 0, dim - 1, 1}];
  Subscript[L, k] = 
   Table[Table[L[i, j, k], {j, 0, dim - 1, 1}], {i, 0, dim - 1, 1}];
  Subscript[Ld, k] = 
   Table[Table[Conjugate[L[j, i, k]], {j, 0, dim - 1, 1}], {i, 0, 
     dim - 1, 1}];, {k, 1, lnlindblad}];
Do[
  Do[
    \[Rho][i, j] = Subscript[\[Rho]s, i, j][t];
    , {j, 0, dim - 1, 1}];
  , {i, 0, dim - 1, 1}];
\[Rho] = 
  Table[Table[\[Rho][i, j], {j, 0, dim - 1, 1}], {i, 0, dim - 1, 1}];
Do[
  Do[
    H[i, j] = Subscript[Hs, i, j];
    , {j, 0, dim - 1, 1}];
  , {i, 0, dim - 1, 1}];
H = Table[Table[H[i, j], {j, 0, dim - 1, 1}], {i, 0, dim - 1, 1}];
h = 1;
(*H*)
Table[Subscript[Hs, i, j] = HH[[i + 1, j + 1]], {i, 0, dim - 1}, {j, 
   0, dim - 1}];
Do[Subscript[\[Gamma], i] = gs[[i]], {i, 1, lnlindblad}];
Do[Table[
   Subscript[L, i, j, k] = CCs[[k]][[i + 1, j + 1]], {i, 0, 
    dim - 1}, {j, 0, dim - 1}], {k, 1, lnlindblad}];
LIND = Flatten[
   Expand[(-I/h) Conm[H, \[Rho]] + 
     Sum[Subscript[\[Gamma], 
        k] (Subscript[L, k] . \[Rho] . 
          Subscript[Ld, k] - (1/2) AntConm[
           Subscript[Ld, k] . Subscript[L, k], \[Rho]]), {k, 1, 
       lnlindblad}]]];
LINDBLAD = 
  1. Chop[Normal@CoefficientArrays[LIND, Flatten[\[Rho]]]][[2]];
re = Re[LINDBLAD];
im = Im[LINDBLAD];
LL = Table[{Flatten[re, 1][[i]], Flatten[im, 1][[i]]}, {i, 1, 
    Length[Flatten[re, 1]]}];
Export["L.dat", LL];
Print["LINBLAD Operator crated L.dat"]
Exit[];
