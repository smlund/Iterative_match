(* 
im_seed.m    

Define functions needed to calculate the seed iteration.  

Steve Lund and Sven Chilton
Lawrence Berkeley National Lab 
Jan, 2006 

Updated by Steve Lund and Kei Fukushima (Hiroshima University) 
June, 2012  
Updated by Steve Lund 
March, 2014; Dec 2016   

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
Seed[]
Function to calculate seed iteration as a function of:
  
  SolCase  .... solution case: 1, and 2 implemented 
  Q        .... perveance [1]
  emitx    .... x rms edge emittance [m-rad]
  emity    .... y rms edge emittance [m-rad]
  sigmax   .... x plane depressed phase advance [rad/period]
  sigmay   .... y plane depressed phase advance [rad/period]

It is assumed that the lattice and continuous limit quantities have 
been appropriately defined (input as global variables) before 
Seed[] is called. The function returns as global variables:

  Q0      .... Q^0 seed dimensionless perveance [1]
  emitx0  .... emit_x^0 seed x-emittance [m-rad]
  emity0  .... emit_y^0 seed y-emittance [m-rad]
  sigmax0 .... sigma_x^0 seed x-plane depressed emittance [rad/period]
  sigmay0 .... sigma_y^0 seed y-plane depressed emittance [rad/period]

  betax0[s] .... beta_x^0(s) seed x-plane beta func on s = [0,lperiod] [m]
  betay0[s] .... beta_y^0(s) seed y-plane beta func on s = [0,lperiod] [m]
 
  rx0[s]  .... r_x^0(s) seed x-plane envelope func on s = [0,lperiod] [m]
  ry0[s]  .... r_y^0(s) seed y-plane envelope func on s = [0,lperiod] [m]

*)

Seed[SolCase_,Q_?NumericQ,emitx_?NumericQ,emity_?NumericQ,
     sigmax_?NumericQ,sigmay_?NumericQ
    ] := 
Module[{Qx0,Qy0,emitr0}, 
  Which[
    SolCase === -1,
      (* --- nothing need to set: direct integration from initial condition *)
      True, 
    SolCase === 1,
      (* --- set known parameters *)
      Q0 = Q;
      sigmax0 = sigmax;
      sigmay0 = sigmay;
      (* --- calculate betatron functions *) 
      focusingx[s_]:=kappax[s]-2*Qbar/((rxbar+rybar)*rxbar);
      focusingy[s_]:=kappay[s]-2*Qbar/((rxbar+rybar)*rybar);
      {betax0[s_],betay0[s_]} = Betafunction[sigmax0,sigmay0,s,focusingx,focusingy,1];
      (* --- calculate emittances from known parameters and 
             betatron functions.  Zero perveance is taken as a special case *)
      emitr0 = ( Bar[kappax[s]*Sqrt[betax0[s]]] - Bar[1/betax0[s]^(3/2)] )/
               ( Bar[kappay[s]*Sqrt[betay0[s]]] - Bar[1/betay0[s]^(3/2)] );
      emitx0 = 2*Q0*Bar[1/( Sqrt[betax0[s]] + Sqrt[emitr0]*Sqrt[betay0[s]] )]/
                    ( Bar[kappax[s]*Sqrt[betax0[s]]] - Bar[1/betax0[s]^(3/2)] );
      emity0 = 2*Q0*Bar[1/( Sqrt[betax0[s]]/Sqrt[emitr0] + Sqrt[betay0[s]] )]/
                    ( Bar[kappay[s]*Sqrt[betay0[s]]] - Bar[1/betay0[s]^(3/2)] );
      (* --- check that both emittances are positive.  if not, use continuous 
             estimates *) 
      If[(emitx0 < 0) || (emity0 < 0), 
         emitx0 = Sqrt[ ( (sigma0x/lperiod)^2*rxbar - 2*Q0/(rxbar+rybar) )*
                        rxbar^3
                      ]; 
         emity0 = Sqrt[ ( (sigma0y/lperiod)^2*rybar - 2*Q0/(rxbar+rybar) )*
                        rybar^3
                      ]
        ], 
    SolCase === 2 || SolCase === 0,
      (* --- set known parameters *)
      emitx0 = emitx;
      emity0 = emity;
      sigmax0 = sigmax;
      sigmay0 = sigmay;
      (* --- calculate betatron functions *)
      focusingx[s_]:=kappax[s]-2*Qbar/((rxbar+rybar)*rxbar);
      focusingy[s_]:=kappay[s]-2*Qbar/((rxbar+rybar)*rybar);
      {betax0[s_],betay0[s_]} = Betafunction[sigmax0,sigmay0,s,focusingx,focusingy,1];
      (* --- Calculate Q0 from known parameters and betatron functions.  
             Use the average of the x- and y-planes *)
      Qx0 = (emitx0/2)*
        ( Bar[kappax[s]*Sqrt[betax0[s]]] - Bar[1/betax0[s]^(3/2)] )/ 
          Bar[1/( Sqrt[betax0[s]]+Sqrt[emity0/emitx0]*Sqrt[betay0[s]] )];
      Qy0 = (emity0/2)*
        ( Bar[kappay[s]*Sqrt[betay0[s]]] - Bar[1/betay0[s]^(3/2)] )/ 
          Bar[1/( Sqrt[emitx0/emity0]*Sqrt[betax0[s]]+Sqrt[betay0[s]] )];
      Q0 = ( Qx0 + Qy0 )/2;
      (* --- check that both perveances are positive.  if not, use continuous 
             estimates *) 
      If[(Qx0 < 0) || (Qy0 < 0), 
         Qx0 = ( (sigma0x/lperiod)^2*rxbar - emitx0^2/rxbar^3 )*
               (rxbar+rybar)/2; 
         Qy0 = ( (sigma0y/lperiod)^2*rybar - emity0^2/rybar^3 )*
               (rxbar+rybar)/2
        ], 
    1 === 1, 
      Print["Seed Case Not Defined"];
      Abort[] 
       ];
  (* --- define new iteration *)
  rx0[s_] := Sqrt[emitx0*betax0[s]];
  ry0[s_] := Sqrt[emity0*betay0[s]];
  (* --- dummy return *)
  True
      ];


(* DEBUG:  execute after reading in script to test *) 

If[DeBug,
(* *)
Print["DeBug: im_seed: Seed Iteration Test Results:"];   
ContLimit[SolCase,Q,emitx,emity,sigmax,sigmay];
Seed[SolCase,Q,emitx,emity,sigmax,sigmay];
(* *)
Print["DeBug: im_seed: Matched Depressed Betatron Function Seed Iteration:"];
Print[ 
Plot[ {betax0[s],betay0[s]}, {s,0,lperiod},
  PlotStyle -> {Black,Red},
  PlotLabel -> "Seed Betatron Functions", 
  AxesLabel -> {"s (m)","\!\( \[Beta]\_x \), \!\( \[Beta]\_y \) (m)"}
    ]
     ];   
Print["DeBug: im_seed: Matched Envelope Seed Iteration:"]; 
Print[
Plot[ {rx0[s],ry0[s]}, {s,0,lperiod},
  PlotStyle -> {Black,Red},
  PlotLabel -> "Seed Envelope Functions", 
  AxesLabel -> {"s (m)","\!\( r\_x \), \!\( r\_y \) (m)"}
    ]
     ]
   ];
