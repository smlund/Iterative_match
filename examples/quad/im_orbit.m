(* 
Characteristic orbit plots both depressed and undepressed x-plane orbits 
within a matched beam envelope.  

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

Based on a matchecd beam envelope being 
calculated by im_solver.m  


Caution:
  The orbits appear correct, but orbits do not exactly repeat initial 
conditions for one full betatron oscillation in periodic focusing unless 
sigma is a commencerate fraction of 360 degrees.  At first 
I thought this was a program error, but that is likely a correct answer.  

  sigma_0 = 60 deg sigma/sigma_0 = 0.2 => sigma = 0.2*60=12 
        => orbits repeat 360/12 = 30 periods regardless of how the 
           matched envelope solution (case 0,1,2) is calulated 
  sigma_ 0 = 60 deg, sigma/sigma_0 = 0.6421....   
        => orbits did not exactly repeat after betatron period incommencerate 
           with lattice period
*)

(* Plot Control *)
orbitpltaspectratio = 1/4;   (* -- aspect ratio of plots *)


(* Periodic lattice extensions for matched envelope and focusing funcs 

     Comment:  Had to pad ?NumericQ on arguments to avoid issues when 
     used in NDSolve[] within HillSolverper[].  It appeared that arguments 
     were not being kept in periodic range when passed for matched envelope 
     (for unknown reasons) without this useage.   
*)
rxper[s_?NumericQ] :=
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  rx[ss]
   ];
ryper[s_?NumericQ] := 
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  ry[ss]
   ];
rxpper[s_?NumericQ] := 
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  rxp[ss]
   ];
rypper[s_?NumericQ] := 
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  ryp[ss]
   ];
kappaxper[s_?NumericQ] := 
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  kappax[ss]
   ];
kappayper[s_?NumericQ] := 
Module[{ss}, 
  ss = soi + Mod[s-soi,lperiod];
  kappay[ss]
   ];

(* 
Plot[rx[s], {s,0,lperiod}]
Plot[ rxper[s], {s,0,2*lperiod}]
*)


(* Set initial conditions *)

(* --- initial condition generation for angles consistent with 
       oscillation amp and moving within KV envelope 
*)
kvshellanglex[frac_,x_,rx_,rxp_] := rxp*(x/rx) + 
  (emitx/rx)*Sqrt[Abs[frac^2 - (x/rx)^2]]; 

kvshellangley[frac_,y_,ry_,ryp_] := ryp*(y/ry) + 
  (emity/ry)*Sqrt[Abs[frac^2 - (y/ry)^2]];

If[ xoi  === Fenvx,   xoi = fenvx*rxper[soi]];
If[ xpoi === Femitx, 
      xpoi = Chop[ kvshellanglex[femitx,xoi,rxper[soi],rxpper[soi]] ]
  ]; 

If[ yoi  === Fenvy,   yoi = fenvy*ryper[soi]];
If[ ypoi === Femity,  
      ypoi = Chop[ kvshellangley[femity,yoi,ryper[soi],rypper[soi]] ]
  ]; 

(* --- reset final orbit coordinate if specified for period *)

If[ sof === Period, sof = soi + (2*Pi/sigmax)*lperiod]; 


(* Calculate orbits *)

(* --- extended period Hill's equation solver 
       Note brkptsx and brkptsy are defined on [0,lperiod] 
*)

nperiodint = IntegerPart[nperiod + 0.5];
brkptsperx = Flatten[Table[soi + i*lperiod + brkptsx, {i, 0, nperiodint-1}]];
brkptspery = Flatten[Table[soi + i*lperiod + brkptsy, {i, 0, nperiodint-1}]];

HillSolveper[kappa_,xi_?NumericQ,xpi_?NumericQ,brkptsper_] :=  
Module[{dummy},
  (* --- solve differential equation *)
  NDSolve[
    {xp[s] == x'[s],
     xp'[s] + kappa[s]*x[s] == 0,
     x[soi] == xi, xp[soi] == xpi
    }, 
    {x,xp}, {s,soi,brkptsper,sof},
    MaxStepSize   -> lperiod/200,
    AccuracyGoal  -> ndag, 
    PrecisionGoal -> ndpg
         ]
      ];

(* --- solve for x- and y-plane orbits; both depressed and 
       undepressed 
*)

focxper[s_] := kappaxper[s] - 2*Q/(rxper[s]*(rxper[s]+ryper[s]));
xorbitsolve = 
HillSolveper[focxper,xoi,xpoi,brkptsperx];

xo  = x  /. xorbitsolve[[1]];
xpo = xp /. xorbitsolve[[1]];
          
  
(*
Plot[ xo[s], {s,soi,sof}]
Plot[ kappaxper[s], {s,soi,sof}]
*)


focyper[s_] := kappayper[s] - 2*Q/(ryper[s]*(rxper[s]+ryper[s]));
yorbitsolve = 
HillSolveper[focyper,yoi,ypoi,brkptspery];

yo  = x  /. yorbitsolve[[1]];
ypo = xp /. yorbitsolve[[1]];


xorbitsolve0 = 
HillSolveper[kappaxper,xoi,xpoi,brkptsperx];

xo0  = x  /. xorbitsolve0[[1]];
xpo0 = xp /. xorbitsolve0[[1]];

(*
Plot[ xo0[s], {s,soi,sof}]
Plot[ kappaxper[s], {s,soi,sof}]
*)

yorbitsolve0 = 
HillSolveper[kappayper,yoi,ypoi,brkptspery];

yo0  = x  /. yorbitsolve0[[1]];
ypo0 = xp /. yorbitsolve0[[1]];


(* Courant Snyder invariant *)

(* -- periodic w-functions with space-charge *)
wxper[s_]  := rxper[s]/Sqrt[emitx];
wxpper[s_] := rxpper[s]/Sqrt[emitx];

wyper[s_]  := ryper[s]/Sqrt[emity];
wypper[s_] := rypper[s]/Sqrt[emity];

csinv[x_,xp_,wx_,wxp_] := (x/wx)^2 + (wx*xp-wxp*x)^2;


(*
Plot[csinv[xo[s],xpo[s],wxper[s],wxpper[s]], {s,soi,sof}, 
PlotRange -> All]
*)

(* 
Accuracy check ... this takes a long time 
epsx = 
NIntegrate[ csinv[xo[s],xpo[s],wxper[s],wxpper[s]], 
            {s,soi,sof}
          ]/(sof-soi);
*)

epsx = csinv[xo[soi],xpo[soi],wxper[soi],wxpper[soi]];
epsy = csinv[yo[soi],ypo[soi],wyper[soi],wypper[soi]];

StylePrint["Characteristic x- and y-Plane Orbits","Section"];


(* --- information printouts *)
headingorbits = {
{"Single Particle CS Invariants (includes space-charge):", 
 "  \!\( \[Epsilon]\_x \) [mm-mrad]",
 "  \!\( \[Epsilon]\_y \) [mm-mrad]",
 "Axial Coordinates:", 
 "  Initial \!\( s\_i \) [m]",
 "  Final   \!\( s\_f \) [m]", 
 "Initial Conditions, Undep and Dep",
 "   x-plane",
 "     x[\!\( s\_i \)] [mm]", 
 "     x'[\!\( s\_i \)] [mrad]",
 "   y-plane",
 "     y[\!\( s\_i \)] [mm]", 
 "     y'[\!\( s\_i \)] [mrad]", 
 "Final Conditions, Undepressed", 
 "   x-plane",
 "     x[\!\( s\_f \)] [mm]", 
 "     x'[\!\( s\_f \)] [mrad]",
 "   y-plane",
 "     y[\!\( s\_f \)] [mm]", 
 "     y'[\!\( s\_f \)] [mrad]",
 "Final Conditions, Depressed", 
 "   x-plane",
 "     x[\!\( s\_f \)] [mm]", 
 "     x'[\!\( s\_f \)] [mrad]",
 "   y-plane",
 "     y[\!\( s\_f \)] [mm]", 
 "     y'[\!\( s\_f \)] [mrad]" 
},None         };

outputorbits = {
 " ",
 NumberForm[N[1000^2*epsx], ndigits],
 NumberForm[N[1000^2*epsy], ndigits],
 " ", 
 NumberForm[N[soi], ndigits],
 NumberForm[N[sof], ndigits], 
 " ", 
 " ", 
 NumberForm[1000*xo[soi], ndigits], 
 NumberForm[1000*xpo[soi],ndigits],
 " ", 
 NumberForm[1000*yo[soi], ndigits], 
 NumberForm[1000*ypo[soi],ndigits], 
 " ", 
 " ", 
 NumberForm[1000*xo0[sof], ndigits], 
 NumberForm[1000*xpo0[sof],ndigits],
 " ", 
 NumberForm[1000*yo0[sof], ndigits], 
 NumberForm[1000*ypo0[sof],ndigits], 
 " ", 
 " ", 
 NumberForm[1000*xo[sof], ndigits], 
 NumberForm[1000*xpo[sof],ndigits],
 " ", 
 NumberForm[1000*yo[sof], ndigits], 
 NumberForm[1000*ypo[sof],ndigits]
                };

StylePrint[ TableForm[outputorbits,TableHeadings -> headingorbits] ];




(*  
Make x- and y-plane orbit plots with plus/minus matched envelope 
plane projection. 
*)

xorbitplt = 
Plot[ {1000*rxper[ss*lperiod],-1000*rxper[ss*lperiod],
       1000*xo[ss*lperiod],1000*xo0[ss*lperiod]}, 
      {ss,soi/lperiod,sof/lperiod},
  PlotRange       -> {{soi/lperiod,sof/lperiod},All},
  PlotStyle       -> {Blue,Blue,Black,Red},
  Axes            -> False, 
  Frame           -> True, 
  PlotPoints      -> IntegerPart[200*nperiod], 
  AspectRatio     -> orbitpltaspectratio, 
  FrameLabel      -> {"s/\!\( L\_p \), Lattice Periods","x [mm]",
                      "x-Orbits: Red=Undep, Black=Dep, Blue = pm Env",None},
  ImageSize       -> plotwidth,
  BaseStyle       -> basestyle
    ];

yorbitplt = 
Plot[ {1000*ryper[ss*lperiod],-1000*ryper[ss*lperiod],
       1000*yo[ss*lperiod],1000*yo0[ss*lperiod]}, 
      {ss,soi/lperiod,sof/lperiod},
  PlotRange       -> {{soi/lperiod,sof/lperiod},All},
  PlotStyle       -> {Blue,Blue,Black,Red},
  Axes            -> False, 
  Frame           -> True, 
  PlotPoints      -> IntegerPart[200*nperiod], 
  AspectRatio     -> orbitpltaspectratio, 
  FrameLabel      -> {"s/\!\( L\_p \), Lattice Periods","y [mm]",
                      "y-Orbits: Red=Undep, Black=Dep, Blue = pm Env",None},
  ImageSize       -> plotwidth,
  BaseStyle       -> basestyle
    ];

Print[xorbitplt];
Print[yorbitplt];
 
Export["im_plt_orbit_x.eps",xorbitplt,"EPS"];
Export["im_plt_orbit_y.eps",yorbitplt,"EPS"];


(*  
Make x' and y' orbit plots. 
*)

xporbitplt = 
Plot[ {1000*xpo[ss*lperiod],1000*xpo0[ss*lperiod]}, 
      {ss,soi/lperiod,sof/lperiod},
  PlotRange       -> {{soi/lperiod,sof/lperiod},All},
  PlotStyle       -> {Black,Red},
  Axes            -> False, 
  Frame           -> True, 
  PlotPoints      -> IntegerPart[200*nperiod], 
  AspectRatio     -> orbitpltaspectratio, 
  FrameLabel      -> {"s/\!\( L\_p \), Lattice Periods","x' [mrad]",
                      "x' Orbits: Red=Undep, Black=Dep",None},
  ImageSize       -> plotwidth,
  BaseStyle       -> basestyle
    ];

yporbitplt = 
Plot[ {1000*ypo[ss*lperiod],1000*ypo0[ss*lperiod]}, 
      {ss,soi/lperiod,sof/lperiod},
  PlotRange       -> {{soi/lperiod,sof/lperiod},All},
  PlotStyle       -> {Black,Red},
  Axes            -> False, 
  Frame           -> True, 
  PlotPoints      -> IntegerPart[200*nperiod], 
  AspectRatio     -> orbitpltaspectratio, 
  FrameLabel      -> {"s/\!\( L\_p \), Lattice Periods","y' [mrad]",
                      "y' Orbits: Red=Undep, Black=Dep",None},
  ImageSize       -> plotwidth,
  BaseStyle       -> basestyle
    ];

Print[xporbitplt];
Print[yporbitplt];
 
Export["im_plt_orbit_xp.eps",xporbitplt,"EPS"];
Export["im_plt_orbit_yp.eps",yporbitplt,"EPS"];
