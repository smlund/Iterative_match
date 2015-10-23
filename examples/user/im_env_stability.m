(* 
Matched envelope linear stability properties    

Based on a matchecd beam envelope being calculated by im_solver.m.  Stability 
properties are reported within the framework discussed in the paper:

Lund and Bukh, "Stability properties of the tranverse envelope equations 
describing intense ion beam transport," Phys. Rev. Special Topics -- 
Accelerators and Beams 7, 024801 (2004)

Updated by Steve Lund and Kei Fukushima (Hiroshima University) 
June, 2012
Updated by Steve Lund 
March, 2014  

Contact:
Steven M. Lund
Facility for Rare Isotope Beams 
Michigan State University 
640 South Shaw Lane 
East Lansing, MI 48824 
 
lund@frib.msu.edu 
517-908-7291 

*)

(* 
Numerical control parameters 
*)
fuzz = 10^-5;

(* 
Continuous limit mode phase advances
  sigmap = breathing   (plus)  mode 
  sigmam = quadrupoole (minus) mode  

  Average between planes for accuracy in symmetric case and for an approx 
  definition in more general lattices 
*)
sigmapx = Sqrt[2*sigma0x^2 + 2*sigmax^2];
sigmamx = Sqrt[sigma0x^2 + 3*sigmax^2];

sigmapy = Sqrt[2*sigma0y^2 + 2*sigmay^2];
sigmamy = Sqrt[sigma0y^2 + 3*sigmay^2];

sigmap = (sigmapx+sigmapy)/2;
sigmam = (sigmamx+sigmamy)/2;


(* 
Calculate homogeneous transfer matrix of envelope perturbaitions through 
one lattice period 
*)

(* --- ODE advance of envelope perturbations from an initial condition 
       through one lattice period 
*)

pertenvsol[drxi_,drxpi_,dryi_,drypi_] := 
NDSolve[ 
  {drxp'[s] + kappax[s]*drx[s] + (2*Q/(rx[s]+ry[s])^2)*(drx[s]+dry[s]) + 
     (3*emitx^2/(rx[s])^4)*drx[s] == 0, 
   drx'[s] == drxp[s], 
   dryp'[s] + kappay[s]*dry[s] + (2*Q/(rx[s]+ry[s])^2)*(drx[s]+dry[s]) + 
     (3*emity^2/(ry[s])^4)*dry[s] == 0, 
   dry'[s] == dryp[s], 
   drx[0] == drxi, drxp[0] == drxpi, 
   dry[0] == dryi, dryp[0] == drypi
  }, 
  {drx,drxp,dry,dryp}, {s,0,lperiod} 
       ];


(* --- Use perturbation solutions to form transfer matrix for 
       dR = (drx,drx',dry,dry') 
*)

EnvPertTransMatrix = 
Module[{m11,m12,m13,m14,
        m21,m22,m23,m24,
        m31,m32,m33,m34,
        m41,m42,m43,m44,
        row1,row2,row3,row4
       },
  (* column solutions *)
  (* -- 1 *)
  row1 = pertenvsol[1,0,0,0];
  m11 = drx[lperiod]  /. row1[[1]];
  m21 = drxp[lperiod] /. row1[[1]];
  m31 = dry[lperiod]  /. row1[[1]];
  m41 = dryp[lperiod] /. row1[[1]];
  (* -- 2 *)
  row2 = pertenvsol[0,1,0,0];
  m12 = drx[lperiod]  /. row2[[1]];
  m22 = drxp[lperiod] /. row2[[1]];
  m32 = dry[lperiod]  /. row2[[1]];
  m42 = dryp[lperiod] /. row2[[1]];
  (* -- 3 *)
  row3 = pertenvsol[0,0,1,0];
  m13 = drx[lperiod]  /. row3[[1]];
  m23 = drxp[lperiod] /. row3[[1]];
  m33 = dry[lperiod]  /. row3[[1]];
  m43 = dryp[lperiod] /. row3[[1]];
  (* -- 4 *)
  row4 = pertenvsol[0,0,0,1];
  m14 = drx[lperiod]  /. row4[[1]];
  m24 = drxp[lperiod] /. row4[[1]];
  m34 = dry[lperiod]  /. row4[[1]];
  m44 = dryp[lperiod] /. row4[[1]];
  (* return constructed transfer matrix *) 
       {{m11,m12,m13,m14}, 
        {m21,m22,m23,m24}, 
        {m31,m32,m33,m34}, 
        {m41,m42,m43,m44}       }   
      ];

(* Calculate eigenvalues and eigenvectors of transfer matrix *)

{evalues,evectors} = Eigensystem[EnvPertTransMatrix];

(* --- list of eigenvalue properties 
       evalmag = magnitudes 
       evalarg = argument (angles) in radians: -Pi to Pi 
       evalang = angles in radians: 0 to 2*Pi
*)

AngleMap[x_] := If[x < 0, 2*Pi + x, x];


evalmag = Abs[evalues];   
evalarg = Arg[evalues];
evalang = Map[AngleMap,evalarg];



(* 
Determine envelope mode symmetry classes based on Lund and Bukh paper. 
Use fuzz factor to account for finite numerical precision.   
*)

Which[ 
  (* --- all on unit circle => Class A *)
  Table[ Abs[evalmag[[i]]-1] < fuzz, {i,1,4}] == {True,True,True,True}, 
    sym = "Class A", 
  (* --- All on real axis => Class D *)
  Table[ Abs[Im[evalues[[i]]]] < fuzz, {i,1,4}] == {True,True,True,True}, 
    sym = "Class D", 
  (* --- All off unit circle and not Class D, then Class B *)
  Table[ Abs[evalmag[[i]]-1] < fuzz, {i,1,4}] ==  {False,False,False,False},
    sym = "Class B", 
  (* --- If not Class A,B,D, then must be Class C *)
  1 == 1, 
    sym = "Class C"
     ];

(*  Info on mode symmetries *)
 

If[sym == "Class A", 
 modeangles = Sort[evalang,Greater];
 headermodes = 
{
"  Mode 1:", 
"    \!\(\[Sigma]\_1 \) [deg/period]", 
"    \!\(\[Gamma]\_1 \)", 
"  Mode 2:",
"    \!\(\[Sigma]\_2 \) [deg/period]", 
"    \!\(\[Gamma]\_2 \)"
};
outputmodes = 
{
" ", 
NumberForm[modeangles*180/Pi,ndigits],
"1",
" ",
NumberForm[modeangles*180/Pi,ndigits],
"1"
}
  ];
If[sym == "Class B", 
 modea = Sort[evalang,Greater];
 modeangles = {modea[[1]],modea[[4]]};
 modegam    = Sort[evalmag,Greater][[1]];
 headermodes = 
{
"  Mode 1:", 
"    \!\(\[Sigma]\_1 \) [deg/period]", 
"    \!\(\[Gamma]\_1 \)", 
"  Mode 2:",
"    \!\(\[Sigma]\_2 \) [deg/period]", 
"    \!\(\[Gamma]\_2  = 1/\[Gamma]\_1 \)"
};
outputmodes = 
{
" ", 
NumberForm[modeangles*180/Pi,ndigits],
NumberForm[modegam,ndigits],
" ",
NumberForm[modeangles*180/Pi,ndigits],
NumberForm[1/modegam,ndigits]
}
  ];
If[sym == "Class C", 
 evalchop = Chop[evalues];
 evalreal = Select[evalchop,Head[#] === Real &];
 evalreal = Sort[evalreal,Greater];
 gamma1 = evalreal[[1]];
 evalcomp = Select[evalchop,Head[#] === Complex &];
 evalcomparg = Map[AngleMap,Arg[evalcomp]];
 headermodes = 
{
"  Mode 1:", 
"    \!\(\[Sigma]\_1 \) [deg/period]", 
"    \!\(\[Gamma]\_1 \)", 
"  Mode 2:",
"    \!\(\[Sigma]\_2 \) [deg/period]", 
"    \!\(\[Gamma]\_2 \)",
"  Mode 3:", 
"    \!\(\[Sigma]\_3 \) [deg/period]", 
"    \!\(\[Gamma]\_3 = 1/\[Gamma]\_2\)"
};
outputmodes = 
{
" ", 
NumberForm[evalcomparg*180/Pi,ndigits],
"1",
" ",
"180",
NumberForm[gamma1,ndigits],
" ", 
"180",
NumberForm[1/gamma1,ndigits]
}
  ]; 
If[sym == "Class D",
 evalmagsort = Sort[evalmag,Greater];
 gamma1 = evalmagsort[[1]];
 gamma2 = evalmagsort[[2]];
 headermodes = 
{
"  Mode 1:", 
"    \!\(\[Sigma]\_1 \) [deg/period]", 
"    \!\(\[Gamma]\_1 \)", 
"  Mode 2:",
"    \!\(\[Sigma]\_2 \) [deg/period]", 
"    \!\(\[Gamma]\_2 \)",
"  Mode 3:", 
"    \!\(\[Sigma]\_3 \) [deg/period]", 
"    \!\(\[Gamma]\_3 = 1/\[Gamma]\_1\)",
"  Mode 4:",  
"    \!\(\[Sigma]\_4 \) [deg/period]", 
"    \!\(\[Gamma]\_4 = 1/\[Gamma]\_2\)"
};
outputmodes = 
{
" ", 
"180",
NumberForm[gamma1,ndigits],
" ",
"180",
NumberForm[gamma2,ndigits],
" ", 
"180",
NumberForm[1/gamma1,ndigits],
" ",
"180",
NumberForm[1/gamma2,ndigits]
}
  ];


(* Output results *)

StylePrint["Envelope Linear Stability","Section"];


(* --- information printouts *)
headingenvstab = {
{"Continuous Limit Mode Phase Advances:",
 "  (x-y plane averages)", 
 "  \!\( \[Sigma]\_+ \) [deg/period]",
 "  \!\( \[Sigma]\_- \) [deg/period]", 
 "Linear Eigenvalues {|\!\(\[Lambda]\)|,Arg[\!\(\[Lambda]\)]}  {[1],[deg]}:", 
 "  \!\(\[Lambda]\_1\)",
 "  \!\(\[Lambda]\_2\)",
 "  \!\(\[Lambda]\_3\)",
 "  \!\(\[Lambda]\_4\)", 
 "Mode Symmetry [Lund and Bukh, PRSTAB (2004)]: "<>sym, 
 "Eigen Modes:" 
},None         };

outputenvstab = { 
 " ", 
 " ", 
 NumberForm[N[(180/Pi)*sigmap], ndigits],
 NumberForm[N[(180/Pi)*sigmam], ndigits], 
 " ", 
 {NumberForm[evalmag[[1]],ndigits],NumberForm[(180/Pi)*evalarg[[1]],ndigits]},
 {NumberForm[evalmag[[2]],ndigits],NumberForm[(180/Pi)*evalarg[[2]],ndigits]},
 {NumberForm[evalmag[[3]],ndigits],NumberForm[(180/Pi)*evalarg[[3]],ndigits]},
 {NumberForm[evalmag[[4]],ndigits],NumberForm[(180/Pi)*evalarg[[4]],ndigits]},
 " ", 
 " "
                };

StylePrint[ TableForm[outputenvstab,TableHeadings -> headingenvstab] ];
StylePrint[ TableForm[outputmodes,TableHeadings -> {headermodes,None}] ];


(* --- envelope eigenvalue plots on complex unit circle *)

ucirc = Graphics[ Circle[{0,0},1] ];
maxextent = 1.5
axes  = Graphics[ Line[{{{-maxextent,0},{maxextent,0}},
                        {{0,-maxextent},{0,maxextent}}
                        }
                      ] 
                ];
evpts = Graphics[ 
 {PointSize[Large], 
  Point[Table[{Re[evalues[[i]]],Im[evalues[[i]]]}, {i,1,4}]]
 } 
                ];
envevalplt = 
Show[ucirc,axes,evpts,
  ImageSize   -> plotwidth/2, 
  PlotLabel   -> "Linear Perturbation Eigenvalues", 
  BaseStyle   -> basestyle
    ];

Print[envevalplt];

Export["im_plt_env_stab_evals.eps",envevalplt,"EPS"];

