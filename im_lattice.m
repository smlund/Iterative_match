(* 
im_lattice.m

Define piecewise constant lattice focusing functions with the correct 
scale needed to achieve the specified undepressed phase advance values.  

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
General comments:

The user can define their own lattice focusing functions:

   kappax[s]  .... kappa_x(s) 
   kappay[s]  .... kappa_y(s) 

in this module to work with their specific lattices of interest.   
In making these definitions, the functions should be defined as 
numerical in Mathematica for optimal program efficiency using 

   kappax[s_?NumericQ] := ..... 

The amplitude of the kappax[s] and kappay[s] functions should be set 
consistently with the undepressed phase advances sigma0x and sigma0y.  
Or conversely, if specific amplitudes of the kappax[s] and kappay[s] 
functions are employed, the undepressed phase advances should be 
set consistently.  This can be done by analyzing single particle 
transfer matrices of a particle in the applied field advanced 
through one lattice period as explained in any 
standard accelerator physics books/reviews.  See, for example, 
   
   H. Wiedemann, Particle Accelerator Physics: Basic Principles and 
   Linear Dynamics (Springer-Verlag, New York, 1993).  

   S.M. Lund and B. Bukh, "Stability properties of the transverse envelope 
   equations describing intense ion beam transport," 
   Phys. Rev. Special Topics -- Accelerators and Beams 7, 024801 (2004). 

The max focusing amplitude kappa should also be identified, if desired 
(it is not necessary to do so), for informational output used at the end 
of the matching program.  

In the example focusing functions for continuous, solenoidal, and 
quadrupole doublet lattices implemented below, constraint equations 
derived in Lund and Bukh, PRSTAB 7, 024801 (2004) are employed to set 
the overall scale consistently with specified undepressed phase advance 
parameters.  These constraints are derived based on single particle 
orbit transfer matrices in these simple focusing lattices.  Scale factors  
needed are calculated to high precision to ensure that values obtained 
do not influence the precision of matched envelope solutions subsequently 
calculated.  
*)

(*
Read in functions to set up lattice.  This script will define 
lattice focusing functions and parameters associated with the 
classes of piecewise constant kappa functions considered.   


  d  ....  solenoid drift [meters]
  d1 ....  doublet quadrupole drift 1 [meters] 
  d2 ....  doublet quadrupole drift 2 [meters]

      Note: for LatticeType = UserInput these values are ignored. 

  kappax[s] ... lattice x-plane focusing function [meters^-2]
  kappay[s] ... lattice y-plane focusing function [meters^-2]

  brkptsx .... list of kappax discontinuities in lattice period 
  brkptsy .... list of kappay discontinuities in lattice period 

  s = [0,lperiod] in coordinate system setup to define the lattice 
*)

(*
Piecewise constant solenoid lattice drifts.  
  Used only for LatticeType = Solenoid 
  in setting up kappax and kappay 
*)

d  = (1-eta)*lperiod;

(* 
Piecewise constant quadrupole lattice drifts. 
  Used only for LatticeType = Quadrupole 
  in setting up kappax and kappay 
*)

d1 = alpha*(1-eta)*lperiod;
d2 = (1-alpha)*(1-eta)*lperiod;

(*
Reset lattice value of eta for continuous focusing if not set correctly  
*)

If[ LatticeType === Continuous && (eta != 1), eta = 1];

(*
Define common-plane phase advance for the classes of lattices 
considered 
*)

sigma0 = sigma0x;
If[ sigma0x != sigma0y, 
    Print["Error: sigmax not equal sigmay"]; 
    Abort[] 
  ];

(*
Calculate lattice focusing strength Theta from sigma0  
and numerical root finding of the relevant phase advance 
equation of constraint for the type of lattice specified.  
*)

(* --- Phase advance equation of constraint based on lattice type. *)
cossigma0[theta_] := 
 Which[
   LatticeType === Continuous,
     Cos[2*theta],
   LatticeType === Solenoid,
     Cos[2*theta] - ((1-eta)/eta)*theta*Sin[2*theta],
   LatticeType === Quadrupole,
     Cos[theta]*Cosh[theta] +
       ((1-eta)/eta)*theta*(Cos[theta]*Sinh[theta]-Sin[theta]*Cosh[theta]) -
       (2*alpha)*(1-alpha)*((1-eta)*theta/eta)^2*Sin[theta]*Sinh[theta], 
   1 === 1, 
     Null 
      ];

(* --- Guessed value of theta for use in numerical root finding. 
       The value employed is exact for continuous focusing and is 
       based on thin lens estimates for solenoidal and quadrupole focusing 
*)
thetaguess = 
  Which[
    LatticeType === Continuous, 
      sigma0/2, 
    LatticeType === Solenoid,
      Sqrt[eta*(1-Cos[sigma0])/2],
    LatticeType === Quadrupole,
      (2*(1-Cos[sigma0])*eta^2/(1-2*eta/3))^(1/4), 
    1 === 1, 
      Null 
       ];

(* --- Calculate needed value of Theta for specified phase advance using 
       numerical root finding for classes of peicewise constant lattices 
       defined. Use high fractional accuracy to ensure 
       the value will be more than accurate enough for any reasonable 
       matched envelope solution tolerance specified.  
*)
If[ LatticeType === Continuous || LatticeType === Solenoid || 
    LatticeType === Quadrupole, 
       Theta = 
         FindRoot[cossigma0[theta] == Cos[sigma0], {theta,thetaguess},
                  WorkingPrecision -> 3*wp,
                  AccuracyGoal     -> Infinity,
                  PrecisionGoal    -> 2*wp
                 ][[1]][[2]],
        (* else *) 
        Theta = Null
   ];

(* --- Calculate focusing amplitude for classes of piecewise constant focusing 
       lattice defined.  
*)
If[ LatticeType === Continuous || LatticeType === Solenoid || 
    LatticeType === Quadrupole, 
      kappa = (2*Theta/(eta*lperiod))^2, 
    (* else *) 
      kappa = Null 
  ];

(*
Define lattice focusing functions on the coordinate interval s = [0,lperiod]. 
  kappax[s] = x-plane focus function 
  kappay[s] = y-plane focus function 
*)

kappax[s_?NumericQ] := 
  Which[
    LatticeType === Continuous, 
      kappa, 
    LatticeType === Solenoid,
      Piecewise[{{kappa, d/2 <= s <= (d/2)+eta*lperiod}}, 0],
    LatticeType === Quadrupole,
      Piecewise[{{kappa,  d2/2 <= s <= (d2+eta*lperiod)/2},
                 {-kappa, d1+(d2+eta*lperiod)/2 <= s <= d1+eta*lperiod+d2/2}
                },0
               ],
    LatticeType === UserInput, 
      kappaxinput[s], 
    1 === 1, 
      Print["Error: im_lattice, LatticeType setting not valid"]; 
      Abort[]     
       ];

kappay[s_?NumericQ] := 
  Which[
    LatticeType === Continuous, 
      kappax[s], 
    LatticeType === Solenoid, 
      kappax[s],
    LatticeType === Quadrupole, 
      -kappax[s], 
    LatticeType === UserInput, 
      kappayinput[s], 
    1 === 1, 
      Print["Error: im_lattice, LatticeType setting not valid"]; 
      Abort[]     
       ];


(* 
Define break points in piecewise constant lattice focusing functions 
that are used in various numerical integrations.  
*)
{brkptsx,brkptsy} = 
  Which[
    LatticeType === Continuous,
      {{lperiod/2},{lperiod/2}},   (* break point in middle just for clean code later *) 
    LatticeType === Solenoid,
      {{d/2, d/2+eta*lperiod},{d/2, d/2+eta*lperiod}},
    LatticeType === Quadrupole,
      {{d2/2, (d2+eta*lperiod)/2, d1+(d2+eta*lperiod)/2, d1+eta*lperiod+d2/2},
       {d2/2, (d2+eta*lperiod)/2, d1+(d2+eta*lperiod)/2, d1+eta*lperiod+d2/2}}, 
    LatticeType === UserInput, 
      {brkptsinputx,brkptsinputy}, 
    1 === 1, 
      Print["Error: im_lattice, LatticeType setting not valid"]; 
      Abort[] 
       ];

brkpts = Union[brkptsx,brkptsy]

(* 
Calculate sigma0x and sigm0y for user input lattice functions 
and in every case calculate the undepressed lattice betatron functions:


  
*)

If[
LatticeType === UserInput,
 (* --- Replace values of sigma0x and sigma0y and calculate 
        undepressed betatron functions 
 *)
 {betax[s_],betay[s_]}= Betafunction[sigma0x,sigma0y,s,kappax,kappay,0];
 (* --- find max value of focusing function in period *)
 kappa = Max[FindMaxMin[kappax][[1]],FindMaxMin[kappay][[1]]],
(* else *)
 (* --- calculate undepressed betatron functions *)
 {betax[s_],betay[s_]}= Betafunction[sigma0x,sigma0y,s,kappax,kappay,1];
 ];


(* DEBUG:  execute after reading in script to test *)

If[DeBug,
Print["DeBug: im_lattice: Diagnostic Lattice Function Plot:"];
Print[
Plot[ {kappax[s],kappay[s]}, {s,0,lperiod},
  PlotStyle -> {Black,Red}, 
  PlotLabel -> "Lattice Focusing Functions", 
  AxesLabel -> {"s (m)","\!\( \[Kappa]\_x \), \!\( \[Kappa]\_y \) (1/m^2)"} 
    ]
     ];

Print["DeBug: im_lattice: Diagnostic Undepressed Betatron Function Plot:"];
Print[
Plot[ {betax[s],betay[s]}, {s,0,lperiod},
  PlotStyle -> {Black,Red}, 
  PlotLabel -> "Lattice Focusing Functions", 
  AxesLabel -> {"s (m)","\!\( \[Beta]\_x \), \!\( \[Beta]\_y \)"} 
    ]
     ];
 ];     
  

(* Calculate depressed phase advances sigmax, sigmay from 
   TuneDepx, TuneDepy and sigma0x, sigma0y 
*)

sigmax = TuneDepx*sigma0x;
sigmay = TuneDepy*sigma0y;