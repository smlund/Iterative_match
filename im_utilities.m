(*
im_utilities.m
 
Utilities used in the IM method.  The lattice should be setup before functions 
contained are defined.     
 
Steve Lund and Sven Chilton
Lawrence Berkeley National Lab
Jan, 2006

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
Bar[] 
Period average function taken over the lattice period 

Note that breakpoints are put in the integration range since 
averages will sometimes be taken of quantities with discontinuities 
and this will help ensure accuracy.  Also, quantities averaged 
with bar must be input using s as the function argument, e.g., 
Bar[rx[s]], etc.    
*)
Bar[f_] :=  
  NIntegrate[f, Evaluate[Flatten[{s,si,brkpts,si+lperiod},1]],
    Method        -> nimethod, 
    AccuracyGoal  -> niag,
    PrecisionGoal -> nipg,
    MaxRecursion  -> 15
            ]/lperiod;

(*
HillSolve[] 
Function to solve Hill's equation with input kappa[s] focusing function

  x''(s) + kappa(s)*x(s) = 0 
  x(s=s_i)  = x_i    initial coordinate 
  x'(s=s_i) = x'_i   initial angle 

Inputs:
   
  kappa[s]  .... focusing function kappa(s) that must be defined from
                   s = si to s = sf  [m^-2]
  si        .... initial axial coordinate s = s_i  
  sf        .... final axial coordinate s = s_f 
  xi        .... initial coordinate x_i [m]
  xpi       .... initial angle x'_i [rad] 
  brkpts    .... breakpoints within interval {si,sf} where kappa[s] may 
                 discontinuously change  

Output:

Interpolating functions for:

  x[s]  .... x(s)  on interval s = [s_i,s_f] [m]
  xp[s] .... x'(s) on interval s = [s_i,s_f] [rad]
  
Useage:

  HillSolve[kappa,si,sf,xi,xpi,brkpts] 
  for a defined function kappa[s] and numerical si,sf,xi,xpi

Comments:
 - Older versions had different argument passing structure which seemed 
   to be fragile under Wolfram Mathematica version changes.  This 
   necessitated code changes.  
 - The breakpoint list passed should not be flagged numerical with ?NumericQ 
   Doing so caused trouble.  
 - Use of breakpoints seems to be an undocumented feature of NDSolve.  
   It may prove fragile.  Logically, may need Flatten[{s,si,brkpts,sf}] 
   at some point for correct syntax.  But seems to work fine without.   
 - Had considerable troubles when numerical precision was set different than 
   default values on AccuracyGoal. This should be kept in mind if tuning 
   values.  
 - Recommend running HillSolve[] on simple test cases to debug if problems 
   arise on future version changes.  
*)
  
HillSolve[kappa_,si_?NumericQ,sf_?NumericQ,xi_?NumericQ,xpi_?NumericQ,brkpts_] :=  
Module[{dummy},
  (* --- solve differential equation *)
  NDSolve[
    {xp[s] == x'[s],
     xp'[s] + kappa[s]*x[s] == 0,
     x[si] == xi, xp[si] == xpi
    }, 
    {x,xp}, {s,si,brkpts,sf},
    MaxSteps      -> 10^6,
    Method        -> ndmethod, 
    AccuracyGoal  -> ndag, 
    PrecisionGoal -> ndpg
      ]
        ];

(*
Debug Check

sol = HillSolve[kappax,0,lperiod,1,0,brkpts]

x = x /. sol[[1]][[1]];

Plot[ x[s], {s,0,lperiod}]

*)

  (*
Betafunction[]
Function to calculate x- and y- plane betatron functions as a function of:

  sigmax     .... x-plane phase advance 
  sigmay     .... y-plane phase advance 
  s          .... axial coordinate
  focusingx  .... x-plane focusing function 
  focusingy  .... y-plane focusing function
  mode       .... 0 : overwrite sigma0x and sigma0y calculated by transfer matrix
                  1 : use input sigmax and sigmay

The function returns {betax[s],betay[s]}.
Here betax[s] and betay[s] are betatron functions of input x-plane and y-plane focusing functions.  

*)

Betafunction[sigmax_?NumericQ,sigmay_?NumericQ,s_,focusingx_,focusingy_,mode_] := 
Module[{solcx1,solsx1,solcy1,solsy1,sigmax1,sigmay1,
        cossigmax1,sinsigmax1,cossigmay1,sinsigmay1,
        cx1,cxp1,sx1,sxp1,cy1,cyp1,sy1,syp1
       },
(* --- Solve differential equations for principal orbits *)
solcx1 = HillSolve[focusingx,si,si+lperiod,cxi,cxpi,brkptsx];
solsx1 = HillSolve[focusingx,si,si+lperiod,sxi,sxpi,brkptsx];
solcy1 = HillSolve[focusingy,si,si+lperiod,cyi,cypi,brkptsy];
solsy1 = HillSolve[focusingy,si,si+lperiod,syi,sypi,brkptsy];
(* --- Extract principal function and derivations into interpolating functions 
       Notation:
         cx1  => c_x^i(s) 
         cxp1 => d c_x^i(s) /ds 
         etc.  
*)
cx1  = x  /.  solcx1[[1]][[1]];
cxp1 = xp /.  solcx1[[1]][[2]];
sx1  = x  /.  solsx1[[1]][[1]];
sxp1 = xp /.  solsx1[[1]][[2]];
cy1  = x  /.  solcy1[[1]][[1]];
cyp1 = xp /.  solcy1[[1]][[2]];
sy1  = x  /.  solsy1[[1]][[1]];
syp1 = xp /.  solsy1[[1]][[2]];
(* --- calculate sin and cos of depressed phase advances *)
Which[
 mode == 0, 
  sigma0x = ArcCos[(cx1[lperiod]+sxp1[lperiod])/2];
  If[Not[Element[sigma0x,Reals]],
     Print["Routine Betafunction: program cannot work with sigma_x complex"];
     Abort[]
    ];
  sigmax1 = sigma0x;
  sigma0y = ArcCos[(cy1[lperiod]+syp1[lperiod])/2];
  If[Not[Element[sigma0y,Reals]],
     Print["Routine Betafunction: program cannot work with sigma_y complex"];
     Abort[]
    ];
  sigmay1 = sigma0y,
 mode == 1, 
  sigmax1 = sigmax; 
  sigmay1 = sigmay
      ];

cossigmax1 = Cos[sigmax1];
sinsigmax1 = Sin[sigmax1];
cossigmay1 = Cos[sigmay1];
sinsigmay1 = Sin[sigmay1];
(* --- calculate and return betatron functions *)

{(* betatron function of x *)
 sinsigmax1*sx1[s]^2/sx1[si+lperiod] + 
               ( sx1[si+lperiod]/sinsigmax1 ) * 
          ( cx1[s] + ((cossigmax1-cx1[si+lperiod])/sx1[si+lperiod])*sx1[s] )^2,
 
 (* betatron function of y *)
 sinsigmay1*sy1[s]^2/sy1[si+lperiod] + 
               ( sy1[si+lperiod]/sinsigmay1 ) * 
          ( cy1[s] + ((cossigmay1-cy1[si+lperiod])/sy1[si+lperiod])*sy1[s] )^2}
      ];


(*
FindMaxMin[]

Find Maximum value and Minimum values with corresponding s-locations 
of an input function defined on the interval s = [0,lperiod].  This works 
by dividing the period into sub-intervals, then running numerical 
optimization on each sub-interval, and returning returning the resulting 
max/min values and locations   

Useage: 

   {funcmax,sfuncmax,funcmin,sfuncmin} = FindMaxMin[func];
 
   func[s] defined on [0,lperiod]

      func[sfuncmax] = funcmax  => max value of function 
      func[sfuncmin] = funcmin  => min vlaue of function 

*)

FindMaxMin[fnc_]:=
Module[{fncper,nsearch,dssearch,fncmax,fncmin,sfncmax,sfncmin,ss,
        findaccuracy,findprecision},

 (* --- define periodically extended function to avoid endpoint issues *)
 fncper[s_] := fnc[Mod[s,lperiod]];

 (* ---- discretized lattice period for subdomains to optimize*)
 nsearch = 50;
 dssearch = lperiod/nsearch;

 (* --- set initial max/mins *) 
 fncmax = fnc[0];
 fncmin = fnc[0];
 sfncmax = 0;
 sfncmin = 0;

 (* --- loop over subdomains to find relative max/min points at subdomain 
        boudndary points
 *) 
 Do[
  ss = dssearch*isearch;
  If[fnc[ss] > fncmax, fncmax = fnc[ss]; sfncmax = ss];
  If[fnc[ss] < fncmin, fncmin = fnc[ss]; sfncmin = ss];
  ,{isearch,1,nsearch,1}
  ];

 (* --- refine max/min values on relevant subdomains using numerical 
        optimization functions. Set accuracy reasonably consistent with 
         
 *)
 findaccuracy = ndigits + 4;
 findprecision = ndigits + 1;

 sfncmax = sx /. 
 FindMinimum[  -fncper[sx], {sx, sfncmax - dssearch/20., sfncmax + dssearch/20.,
                            sfncmax - dssearch, sfncmax + dssearch},
  AccuracyGoal  -> findaccuracy,
  PrecisionGoal -> findprecision
            ][[2]];
 
 sfncmin = sx /. 
 FindMinimum[   fncper[sx], {sx, sfncmin - dssearch/20., sfncmin + dssearch/20.,
                            sfncmin - dssearch, sfncmin + dssearch},
  AccuracyGoal  -> findaccuracy,
  PrecisionGoal -> findprecision
            ][[2]];

  sfncmax = Mod[sfncmax,lperiod]; 
  sfncmin = Mod[sfncmin,lperiod]; 

  fncmax = fnc[sfncmax];
  fncmin = fnc[sfncmin];

  (* --- return results for functions max/min and corresponding s-locations *)
  {fncmax,sfncmax,fncmin,sfncmin}

];



