(*
im_err _chk.m
 
Partial error checks for inputs in IM method program.     
 
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

 
(* LatticeType *)
If[ !(LatticeType === Continuous || LatticeType === Solenoid || 
    LatticeType === Quadrupole || LatticeType === UserInput), 
    Print["Error: LatticeType = " <> ToString[LatticeType] <> " not defined"]; 
    Abort[]
  ];

(* SolCase *)
If[ !(SolCase === 0 || SolCase === 1 || SolCase === 2), 
    Print["Error: SolCase = " <> SolCase <> " not defined"]; 
    Abort[]
  ]; 

(* HybridCase *)
If[ SolCase === 0 && !(HybridCase === 1 || HybridCase === 2), 
    Print["Error: HybridCase = " <> HybridCase <> " not defined"]; 
    Abort[]
  ]; 


(* Undpressed phase advances *) 
(*   --- comment out for now ... unsure if these traps are a problem or not.    
If[ sigma0x > Pi,
    Print["Error: Program not setup to analyze sigma0x > Pi"];
    Abort[] 
  ];
If[ sigma0y > Pi,
    Print["Error: Program not setup to analyze sigma0y > Pi"];
    Abort[] 
  ];
*)

(* Depressed phase advances *)
If[ (SolCase === 1 || SolCase === 2) && sigmax > sigma0x, 
    Print["Error: Parameter input: sigmax > sigma0x with SolCase = 1 or 2"];
    Abort[] 
  ];
If[ (SolCase === 1 || SolCase === 2) && sigmay > sigma0y, 
    Print["Error: Parameter input: sigmay > sigma0y with SolCase = 1 or 2"];
    Abort[] 
  ];
If[ SolCase === 1 && Q == 0 && (TuneDepx != 1.0 || TuneDepy != 1.0), 
    Print["Error: Parameter input: In the case of Q = 0.0,"];
    Print["  tune depressions must be 1.0 with SolCase = 1"];
    Abort[] 
  ];

(* Test for incompatibility of SolCase = 1 and Q = 0 simultaneously. *)
(*
If[ SolCase === 1 && Q == 0.0, 
    Print["Error: Program not setup for SolCase = 1 and Q = 0."];
    Print["  Emittances arbitrary when Q = 0."]; 
    Print["  Try SolCase = 2 with emitx and emity specified"]'
    Abort[]
  ]; 
*)

(* Test for full beam depression in either plane which is incompatible 
   with this version of the IM method *)
If[ SolCase === 0 && (emitx == 0.0 || emity == 0.0), 
    Print["Error: Program not setup for full space-charge depression."];
    Print["  SolCase = 0 with emitx = 0 or emity = 0"];  
    Abort[]
  ]; 
If[ SolCase === 1 && (TuneDepx == 0.0 || TuneDepy == 0.0), 
    Print["Error: Program not setup for full space-charge depression."];
    Print["  SolCase = 1 with sigmax = 0 or sigmay = 0"];  
    Abort[]
  ]; 
If[ SolCase === 2 && (emitx  == 0.0 || emity  == 0.0 || 
                      TuneDepx == 0.0 || TuneDepy == 0.0), 
    Print["Error: Program not setup for full space-charge depression."];
    Print["  SolCase = 2 with one or more of emitx, emity, sigmax, sigmay"];
    Print["  zero."];  
    Abort[]
  ]; 
