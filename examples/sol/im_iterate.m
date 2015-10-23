(* 
im_iterate.m    

Calculate general iteration.  

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
Iterate[] 
Function to calculate updated iteration of betatron and envelope functions 
as a function of:

  SolCase  .... solution case: 1, and 2 implemented 
  Q        .... perveance [1]
  emitx    .... x rms edge emittance [m-rad]
  emity    .... y rms edge emittance [m-rad]
  sigmax   .... x plane depressed phase advance [rad/period]
  sigmay   .... y plane depressed phase advance [rad/period]

It is assumed that the lattice and previous iteration variables  
been appropriately defined (input as global variables) before 
Iterate[] is called by previous calls to Iterate[] or Seed[] 
(initial iteration case).  The function returns as global variables:

  Q1      .... Q^i iteration dimensionless perveance [1]
  emitx1  .... emit_x^i iteration x-emittance [m-rad]
  emity1  .... emit_y^i iteration y-emittance [m-rad]
  sigmax1 .... sigma_x^i iteration x-plane depressed emittance [rad/period]
  sigmay1 .... sigma_y^i iteration y-plane depressed emittance [rad/period]

  betax1[s] .... beta_x^i(s) iteration x-plane beta func on s = [0,lperiod] [m]
  betay1[s] .... beta_y^i(s) iteration y-plane beta func on s = [0,lperiod] [m]
 
  rx1[s]  .... r_x^i(s) iteration x-plane envelope func on s = [0,lperiod] [m]
  ry1[s]  .... r_y^i(s) iteration y-plane envelope func on s = [0,lperiod] [m]

Also returned as a function object:
  
  tol .... fractional difference between updated and previous iteration [1] 

and previous iteration functions are copied to the updated iteration (both 
as global variables) to prepare for further calls to Iterate[]. That is,  
  
  rx0[s] = rx1[s] 
  ry0[s] = ry1[s] 

on exiting the update.  
*)

Iterate[SolCase_,Q_?NumericQ,emitx_?NumericQ,emity_?NumericQ,
        sigmax_?NumericQ,sigmay_?NumericQ
       ] := 
Module[{Qx1,Qy1,emitx1,emity1,emitr1,ittolxlist,ittolylist,ittol}, 
  (* --- solve updated betatron functions *)
  Which[
    SolCase === 1,
      (* --- set known parameters *)
      Q1 = Q;
      Qx1= Q;
      Qy1= Q;
      sigmax1 = sigmax;
      sigmay1 = sigmay;
      (* --- calculate betatron functions *) 
      focusingx[s_]:=kappax[s]-2*Q0/((rx0[s]+ry0[s])*rx0[s]);
      focusingy[s_]:=kappay[s]-2*Q0/((rx0[s]+ry0[s])*ry0[s]);
      {betax1[s_],betay1[s_]} = Betafunction[sigmax1,sigmay1,s,focusingx,focusingy,1];

      (* --- calculate emittances from known parameters and 
             betatron functions *)
      emitr1 = ( Bar[kappax[s]*Sqrt[betax1[s]]] - Bar[1/betax1[s]^(3/2)] )/
               ( Bar[kappay[s]*Sqrt[betay1[s]]] - Bar[1/betay1[s]^(3/2)] );
      emitx1 = 2*Q1*Bar[1/(Sqrt[betax1[s]] + Sqrt[emitr1]*Sqrt[betay1[s]] )]/
                   ( Bar[kappax[s]*Sqrt[betax1[s]]] - Bar[1/betax1[s]^(3/2)] );
      emity1 = 2*Q1*Bar[1/( Sqrt[betax1[s]]/Sqrt[emitr1] + Sqrt[betay1[s]] )]/
                   ( Bar[ kappay[s]*Sqrt[betay1[s]]] - Bar[ 1/betay1[s]^(3/2)] );,
    SolCase === 2 || SolCase === 0,
      (* --- set known parameters *)
      emitx1 = emitx;
      emity1 = emity;
      sigmax1 = sigmax;
      sigmay1 = sigmay;
      (* --- calculate betatron functions *)
      focusingx[s_]:=kappax[s]-2*Q0/((rx0[s]+ry0[s])*rx0[s]);
      focusingy[s_]:=kappay[s]-2*Q0/((rx0[s]+ry0[s])*ry0[s]);
      {betax1[s_],betay1[s_]} = Betafunction[sigmax1,sigmay1,s,focusingx,focusingy,1];
      (* --- Calculate Q0 from known parameters and betatron functions.  
             Use the average of the x- and y-planes *)
      Qx1 = (emitx1/2)*
        ( Bar[kappax[s]*Sqrt[betax1[s]]] - Bar[1/betax1[s]^(3/2)] )/ 
          Bar[1/( Sqrt[betax1[s]]+Sqrt[emity1/emitx1]*Sqrt[betay1[s]] )];
      Qy1 = (emity1/2)*
        ( Bar[kappay[s]*Sqrt[betay1[s]]] - Bar[1/betay1[s]^(3/2)] )/ 
          Bar[1/( Sqrt[emitx1/emity1]*Sqrt[betax1[s]]+Sqrt[betay1[s]] )];
      Q1 = ( Qx1 + Qy1 )/2,
    1 === 1, 
      Print["Seed Case Not Defined"];
      Abort[] 
       ];
  (* --- construct updated iteration envelope functions from 
         betatron functions *)
  rx1[s_] := Sqrt[emitx1*betax1[s]];
  ry1[s_] := Sqrt[emity1*betay1[s]];
  (* ---- sample difference from present to previous iterations to calculate 
          a maximum value of the tolerance achieved *)
  ittolxlist = Abs[(rx1[ssample]-rx0[ssample])/rx1[ssample]];
  ittolylist = Abs[(ry1[ssample]-ry0[ssample])/ry1[ssample]];
  ittol  = Max[ittolxlist,ittolylist];
  (* ---- save envelope iteration and perveance in previous "0" iteration *)
  Q0 = Q1;
  rx0[s_] = rx1[s];
  ry0[s_] = ry1[s];
  (* --- return iteration properties *)
  {ittol,Qx1,Qy1,emitx1,emity1}
      ];


(* 
MatchGen[]
Function to calculate matched solution to a fractional tolerance tol 
that is globally input as a function of the following input parameters:

  SolCase  .... solution case: only cases 1 and 2 are implemented
  SeedGen  .... generate seed iteration to start (True or False) 
                  If False, the seed iteration must already be generated 
                  by some other means such as prior calls to MatchGen[] 
  sigmax   .... x plane depressed phase advance [rad/period]
  sigmay   .... y plane depressed phase advance [rad/period]

Depending on SolCase, global parameters among 

  Q        .... perveance [1]
  emitx    .... x rms edge emittance [m-rad]
  emity    .... y rms edge emittance [m-rad]
 
will also be used in constructing  the matched envelope solution.  The 
parameters not needed under SolCase will be reset to consistent values 
using the applied parameters and the matched envelope solution.   The 
matched envelope functions are returned as global variables:

  rx[s] .... r_x(s) x-plane matched envelope function on s = [0,lperiod] [m]
  ry[s] .... r_x(s) y-plane matched envelope function on s = [0,lperiod] [m]

The function call to MatchGen[] returns:

  {itnum,ittol,Qx1,Qy1,emitx1,emity1} 
    
    itnum  .... number of iterations needed to achieve fractional error tol
    ittol  .... actual tolerance achieved
    Qx1    .... x-plane perveance consistent with final iteration 
    Qy1    .... y-plane perveance consistent with final iteration 
    emitx1 .... x-plane emittance consistent with final iteration 
    emity1 .... y-plane emittance consistent with final iteration 

Also, itnum is used to reset a global variable iternum that is used for 
data structure compatibility issues with other functions.  
*)

MatchGen[SolCase_,SeedGen_,sigmax_?NumericQ,sigmay_?NumericQ] := 
  Module[{ittol=1, itnum=0,itplt0,itplt1},
    (* --- generate seed iteration *)
    If[SeedGen,  
      (* --- generate continuous limit approximation estimates *)
      ContLimit[SolCase,Q,emitx,emity,sigmax,sigmay];
      (* --- generate actual seed iteration *)
      Seed[SolCase,Q,emitx,emity,sigmax,sigmay]
      ];
    (* --- iterate solution till specified tolerance is achieved *)
    While[
      (* --- loop termination criterion *)
      ittol >= tol, 
        (* --- debug info: plot iteration *)
        If[IterationPlot, 
           itplt0 = Plot[{rx0[s],ry0[s]}, {s,0,lperiod}, 
              PlotStyle       -> {{Dashing[{0.005,0.005}],Black},
                                  {Dashing[{0.005,0.005}],Red}
                                 } 
                       ];
          ];
        (* --- calculate trial solution, tolerance, and iterated parameters *)
        {ittol,Qx1,Qy1,emitx1,emity1} = 
          Iterate[SolCase,Q,emitx,emity,sigmax,sigmay];
        (* --- update iteration number *)
        itnum = itnum + 1; 
        (* --- debug info: print iteration data *) 
        If[IterationPrint, 
          Print["iteration ",itnum];
          Print["  tol = ",ScientificForm[ittol]];
          Print["  Qx1 = ",ScientificForm[Qx1]];
          Print["  Qy1 = ",ScientificForm[Qy1]];
          Print["  emitx1 = ",EngineeringForm[emitx1]];
          Print["  emity1 = ",EngineeringForm[emity1]];
          ];
        (* --- debug info: plot iteration *) 
        If[IterationPlot, 
           itplt1 = Plot[{rx1[s],ry1[s]}, {s,0,lperiod}, 
              PlotStyle  -> {Black,Red} 
                        ];
           itplt = 
           Show[itplt1,itplt0,
              PlotLabel  -> "Iteration "<>ToString[itnum]<>
                             " (previous dashed)",
              AxesLabel  -> {"s (m)","r_x, r_y (m)"} 
               ];
           Print[itplt]
          ];
        (* --- terminate iterations if itnum exceeds max value *)
        If[itnum >= itmax, 
          Print["Warning: Max iterations reached without achieving tolerance"];
          Break[] 
          ]
         ];
    (* --- reset any parameters consistent with calculated solution that 
           are NOT (or should not be) specified under the solution case *)
    Which[
      SolCase === 1,
        emitx = emitx1;
        emity = emity1,    
      SolCase === 2, 
        Q = Q1
         ];
    (* --- set matched envelope functions calculated to tolerance *) 
    rx[s_] = rx1[s];
    ry[s_] = ry1[s];
    (* --- set global variable iternum to itnum for compatibility with 
           various function uses *)
    iternum = itnum;
    (* --- return {itnum,ittol,Qx1,Qy1,emitx1,emity1} *)
    {itnum, ittol, Qx1, Qy1, emitx1, emity1} 
   ];

(* 
Match[] 
Function to generate matched envelope solution for all cases specified 
by SolCase.  This works by directly calling MatchGen[] for SolCase = 1 
and SolCase = 2.  For SolCase = 0, numerical root finding combined 
with SolCase = 1 (HybridCase = 1) or SolCase = 2 (HybridCase = 2) type 
iterations to generate the consistent solution.   

Envelope parameters are input as global parameters.  

On exiting, Match[] returns a list:

 {iterations, tolachieved} 

 iterations  .... number iterations needed 
 tolachieved .... actual fractional tolerance achieved
  
Also returned as global variables are the matched envelope functions and 
the envelope parameters.  If a parameter is not required for the 
SolCase specified, then that parameter is ignored in construction of 
the matched solution and then reset to a value consistent with the matched 
envelope (and other parameters used to specify the envelope).

Two auxiliary functions associated with Match[]:
  emitfunc[sigmax,sigmay]    returns {emitx,emity} 
  Qfunc                      returns {Qx,Qy} 
are defined directly below for use in Match[].  These functions are used  
to calculate quantities needed for numerical root finding when implementing 
SolCase = 0 with HybridMethod = 1 or HybridMethod = 2.  These functions 
must be explicitly defined as numerical to be correctly parsed by 
FindRoot.  
*)
emitfunc[sigmax_?NumericQ,sigmay_?NumericQ] := 
  Module[{sol,ex,ey}, 
    sol = MatchGen[1,False,sigmax,sigmay];
    ex = sol[[5]];
    ey = sol[[6]];
    {ex,ey}
        ];
Qfunc[sigmax_?NumericQ,sigmay_?NumericQ] := 
  Module[{sol,Qx,Qy}, 
    sol = MatchGen[2,False,sigmax,sigmay];
    Qx = sol[[3]];
    Qy = sol[[4]];
    {Qx,Qy}
        ];
Match[] := 
  Module[{iterations,tolachieved,sigmaxg,sigmayg,sigroot,itnumlist},
    Which[
      (* --- Q = 0, set rx and ry from the undepressed betafunction 
             functions and emitances.
      *)
      (   SolCase === 0 && Q == 0.0
       || SolCase === 1 && (Q == 0.0 || TuneDepx == 1.0 || TuneDepy == 1.0)
       || SolCase === 2 && (TuneDepx == 1.0 || TuneDepy == 1.0)
      ),
      Q = 0.0;
      sigmax = sigma0x;
      sigmay = sigma0y;
      rx[s_]:=Sqrt[emitx*betax[s]];
      ry[s_]:=Sqrt[emity*betay[s]];
      TuneDepx = 1.;
      TuneDepy = 1.;
      iterations = 0;
      tolachieved = 10^(-ndpg),
      (* --- for cases 1 and 2 directly call MatchGen[] for solution *)
      SolCase === 1 || SolCase === 2, 
        {iterations, tolachieved} = 
          MatchGen[SolCase,True,sigmax,sigmay][[{1,2}]], 
      (* --- for case 0 use hybrid method of numerical root finding combined 
             with case 1 or 2 solutions when Q is not equal to 0.  
      *)
      SolCase === 0 && !(Q == 0), 
        (* --- guess values of sigmax and sigmay from 
               Case 0 continuous limit to seed root finds *)
        ContLimit[SolCase,Q,emitx,emity,sigmax,sigmay];
        sigmaxg = emitx*lperiod/(rxbar^2);
        sigmayg = emity*lperiod/(rybar^2);
        (* --- generate seed for first FindRoot steps *)
        Seed[SolCase,Q,emitx,emity,sigmaxg,sigmayg];
        (* --- initialize list used for iteration counts in FindRoot *) 
        itnumlist = {};
        sigroot = 
          Which[
            (* --- case 1 hybrid solution of case 0 using x,y emittances *)
            HybridCase === 1,
              FindRoot[emitfunc[sx,sy] == {emitx,emity},
                       {{sx, sigmaxg, 1.01*sigmaxg},
                        {sy, sigmayg, 1.01*sigmayg}},
                       Method            -> rfmethod,
                       AccuracyGoal      -> rfag, 
                       PrecisionGoal     -> rfpg, 
                       EvaluationMonitor :> AppendTo[itnumlist,iternum]
                      ],
            (* --- case 2 hybrid solution of case 0 using x,y perveances *)
            HybridCase === 2,
              FindRoot[Qfunc[sx,sy] == {Q,Q},
                       {{sx, sigmaxg, 1.01*sigmaxg},
                        {sy, sigmayg, 1.01*sigmayg}},
                       Method            -> rfmethod,
                       AccuracyGoal      -> rfag, 
                       PrecisionGoal     -> rfpg,
                       EvaluationMonitor :> AppendTo[itnumlist,iternum]
                      ]
               ];
        (* --- extract consistant phase advances from root found solution *)
        sigmax = sx /. sigroot;
        sigmay = sy /. sigroot;
        (* --- count the total number of iterations over all root find steps *)
        iterations = Plus @@ itnumlist;
        (* --- call final iterations with consistent values of sigmax and 
               sigmay to verify that the solution is correct.   This is 
               necessary since the last call to MatchGen[] associated with 
               the FindRoot[] iterations can leave the system in a matched 
               state that is not the desired solution (one with different 
               sigmax and sigmay based on calls last made). We choose 
               SolType = HybridCase for this purpose.  We also do not 
               count these iterations in the total required since they are 
               not strictly necessary if the state of FindRoot[] iterations 
               could be fully controlled.       
        *)
        tolachieved = MatchGen[HybridCase,True,sigmax,sigmay][[2]]
         ];
    (* --- return iterations and achieved tolerance *)
    {iterations, tolachieved} 
        ];

