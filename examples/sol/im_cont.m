(* 
im_cont.m    

Define functions to calculate continuous limit values needed to 
construct a seed iteration.  

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
ContLimit[]
Function that, based on solution type, calculates continuous limit quantities 
as a function of:

  SolCase  .... solution case: 0, 1, and 2 implemented 
  Q        .... perveance [1]
  emitx    .... emit_x x rms edge emittance [m-rad]
  emity    .... emit_y y rms edge emittance [m-rad]
  sigmax   .... sigma_x x plane depressed phase advance [rad/period]
  sigmay   .... sigma_y y plane depressed phase advance [rad/period]

The lattice focusing functions (and sigma0x, sigma0y) are input as 
global variables and must be defined.  The function returns the 
following continuous limit variables as global:
   
  rxbar  .... bar[r_x] x-envelope [m]
  rybar  .... bar[r_y] y-envelope [m]

  Qbar     .... bar[Q] dimensionless perveance [1] 
  emitxbar .... bar[emit_x] x-plane beam rms edge emittance [m-rad]
  emitybar .... bar[emit_y] y-plane beam rms edge emittance [m-rad]

When numerical quantities are calculated, high precision is used to 
ensure that any accuracy limit will not influence subsequent applications.   

Comments:
  In SolCase = 0 
    The structure is general and all limits should work.   
  In SolCase = 1 
    The formulas employed will not work for the zero space-charge limit with 
    sigmaj = sigma0j, or equivalently, Q = 0.   The limit can be analyzed 
    explicitly and the program generalized if this case is desired.  
  In SolCase = 2 
    We assume, for simplicity, that both of the sigmaj and both of the emitj 
    are specified.  This can only be done if the inputs have 
    sigma0x = sigma0y and emitx = emity.  The zero space charge limit should
    work.  
*)

ContLimit[SolCase_,Q_?NumericQ,emitx_?NumericQ,emity_?NumericQ,
          sigmax_?NumericQ,sigmay_?NumericQ
         ] := 
Module[{emitg,sig0g,rguess}, 
  Which[
    SolCase === 0,
      (* --- set known parameters *)
      Qbar = Q;
      emitxbar = emitx;
      emitybar = emity;
      (* --- calculate cont envelope radii with numerical root finding *) 
      emitg = (emitx + emity)/2;
      sig0g = (sigma0x + sigma0y)/2;
      rguess = (lperiod/sig0g)*Sqrt[Qbar/2 + (1/2)*Sqrt[Qbar^2 + 
                4*(sig0g/lperiod)^2*emitg^2]];
      {rxbar, rybar} = {xrx, xry} /. 
         FindRoot[ { (sigma0x/lperiod)^2*xrx - 2*Qbar/(xrx+xry) - 
                       emitxbar^2/xrx^3 == 0,
                     (sigma0y/lperiod)^2*xry - 2*Qbar/(xrx+xry) - 
                     emitybar^2/xry^3 == 0
                   }, {xrx, rguess}, {xry, rguess} (* , Gave Percision trouble
                   WorkingPrecision -> 3*wp,
                   AccuracyGoal     -> Infinity,
                   PrecisionGoal    -> 2*wp        *)
                 ], 
    SolCase === 1,
       (* --- set known parameter *)
       Qbar = Q; 
       (* --- use analytical formulas to calculate cont envelope radii. 
              Treat undepressed case as a special case 
       *)
       rxbar = Sqrt[2*Qbar]*lperiod/Sqrt[(sigma0x^2-sigmax^2) + 
                    (sigma0x^2 - sigmax^2)^2/(sigma0y^2 - sigmay^2)
                     ];

       rybar = Sqrt[2*Qbar]*lperiod/Sqrt[(sigma0y^2-sigmay^2) + 
                    (sigma0y^2 - sigmay^2)^2/(sigma0x^2 - sigmax^2)
                     ];
       (* --- calculate consistent cont focusing emittances *) 
       emitxbar = (sigmax/lperiod)*rxbar^2;
       emitybar = (sigmay/lperiod)*rybar^2,
    SolCase === 2,
       (* --- set known parameters *)
       emitxbar = emitx;
       emitybar = emity;
       (* --- use analytical formulas calculate cont envelope radii *)
       rxbar = Sqrt[emitxbar/(sigmax/lperiod)];
       rybar = Sqrt[emitybar/(sigmay/lperiod)];
       (* --- calculate consistent cont focusing perveance while averaging 
              between two planes to improve accuracy *)
       Qbar = (1/2)*( (sigma0x^2-sigmax^2)*rxbar*(rxbar+rybar)/(2*lperiod^2) + 
                      (sigma0y^2-sigmay^2)*rybar*(rxbar+rybar)/(2*lperiod^2)
                    ) 
       ];
  (* --- dummy return for function *)
  True    
      ]
 

(* DEBUG: execute after reading in script to test  *)

If[DeBug,
Print["DeBug: im_cont: Continuous Lattice Results:"];
ContLimit[SolCase,Q,emitx,emity,sigmax,sigmay];
Print["Qbar  = ",N[Qbar]];
Print["rxbar = ",N[rxbar]," m"];
Print["rybar = ",N[rybar]," m"];
  ];
