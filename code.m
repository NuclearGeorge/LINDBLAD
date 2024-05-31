(* ::Package:: *)
info := {
   Print[
    StyleForm["=====================================================",
      "Section", FontSize -> 14, Black]];
   Print[StyleForm["PACKAGE:", "Section", FontSize -> 14], 
    StyleForm[" LINDBLAD OPERATOR", "Section", FontSize -> 14, 
     Black] ];
   Print[
    StyleForm["BY: Jorge Chávez-Carlos, 2024", "Section", 
     FontSize -> 12, Black]];
   Print[
    StyleForm["=====================================================",
      "Section", FontSize -> 14, Black]];
   Print[
    StyleForm["LINDBLAD_1.2 package:", "Section", FontSize -> 12, 
     Black]];}[[1]];
LINDBLAD[Hi_, Csi_, gsi_, hi_] := {
    Clear[Conm, AntConm, HH, CCs, gss, hh, dim, lnlindblad, 
     CCds, \[Rho], lind, o, LIND, LINDBLADN];
    Conm[A_, B_] := A . B - B . A;
    AntConm[A_, B_] := A . B + B . A;
    HH = Hi;
    CCs = Csi;
    gss = gsi;
    hh = hi;
    dim = Length[HH];
    lnlindblad = Length[CCs];
    CCds = 
     Table[Transpose[Conjugate[CCs[[k]]]], {k, 1, lnlindblad}] ;
    Do[
     Do[
       \[Rho][i, j] = Subscript[\[Rho]s, i, j];
       , {j, 0, dim - 1, 1}];
     , {i, 0, dim - 1, 1}];
    \[Rho] = 
     Table[Table[\[Rho][i, j], {j, 0, dim - 1, 1}], {i, 0, dim - 1, 
       1}];
    LIND = 
     Flatten[Expand[(-I/hh) Conm[HH, \[Rho]] + 
        Sum[gss[[k]] (CCs[[k]] . \[Rho] . 
             CCds[[k]] - (1/2) AntConm[
              CCds[[k]] . CCs[[k]], \[Rho]]), {k, 1, lnlindblad}]]]; 
    LINDBLADN = 
     Simplify[
      1. Chop[Normal@CoefficientArrays[LIND, Flatten[\[Rho]]]][[2]]];
    lind = Flatten[LINDBLADN]; o = OpenWrite["L.dat"];
    Do[Write[o, N[Simplify[lind[[n]]]]], {n, Length[lind]}];
    Close[o];
    Print["Lindbladian operator defined and exported successfully."];
    Clear[Conm, AntConm, HH, CCs, gss, hh, dim, lnlindblad, 
     CCds, \[Rho], lind, o, LIND];
    
    }[[1]];
