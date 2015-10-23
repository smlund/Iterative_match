(* 
im_diag.m

Output information and plots of the matched envelope solution.  

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
Output Properties of Transport Lattice 
*)

(* --- header *)
StylePrint["Transport Lattice","Section"];

(* --- information printouts *)
headinglattice = {
{"Lattice Type",
 "Undepressed Phase Advances [deg/period]",
 "   x-plane, \!\( \[Sigma]\_\(0 x\) \) [deg/period]",
 "   y-plane, \!\( \[Sigma]\_\(0 y\) \) [deg/period]",
 "Lattice Period, \!\( L\_p \) [m]",
 "Occupancy, \!\( \[Eta] \)",
 "Syncopation Factor, \!\(
  \[Alpha] \) (\!\( \[Alpha] = 1/2 \[Implies] \) FODO)",
 "Max Focusing Strength, Max[\!\( \[Kappa]\_x, \[Kappa]\_y \)], [1/\!\( 
  m\^2 \)]"
},None         };

outputlattice = {
 LatticeType,
 " ", 
 NumberForm[N[sigma0x*180/Pi],ndigits],
 NumberForm[N[sigma0y*180/Pi],ndigits],
 NumberForm[N[lperiod],ndigits],
 If[LatticeType===UserInput,"NA",NumberForm[N[eta],ndigits]],
 If[LatticeType===Quadrupole,NumberForm[N[alpha],ndigits],"NA"],
 NumberForm[kappa,ndigits]
                };

StylePrint[ TableForm[outputlattice,TableHeadings -> headinglattice] ];


(* --- make plot of lattice focusing functions *)
PlotLattice = 
Plot[ {kappax[s],kappay[s]}, {s,0,lperiod},
  PlotRange        -> All, 
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> {"s (m)",
   "\!\( \[Kappa]\_x \), \!\( \[Kappa]\_y \) \!\( [m\^\(-2\) ] \)",
                 "Lattice Focusing Functions (black = x, red = y)",None}, 
  ImageSize       -> plotwidth, 
  BaseStyle       -> basestyle   
    ];

(* --- output lattice plot in eps format *)
Export["im_plt_lat_kappa.eps",PlotLattice,"EPS"];


(* --- show lattice plot in notebook *)
Print[PlotLattice];


If[OutputBeta,
(*
Output properties of undepressed betatron function
*)

(* --- header *)
StylePrint["Undepressed (Lattice) Betatron Function","Section"];

(* --- betatron function *)
PlotBeta = 
Plot[ {betax[s],betay[s]}, {s,0,lperiod},
  PlotRange        -> All, 
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> {"s (m)",
   "\!\( \[Beta]\_x \), \!\( \[Beta]\_y \) \!\( [m] \)",
                 "Undepressed Betatron Functions (black = x, red = y)",None}, 
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle  
    ];

    
(* --- output lattice plot in eps format *)
Export["im_plt_lat_beta.eps",PlotBeta,"EPS"];

(* --- show lattice plot in notebook *)
Print[PlotBeta];

(* --- find max/min values of betatron function *)
{betaxmax,sbetaxmax,betaxmin,sbetaxmin}=FindMaxMin[betax];
{betaymax,sbetaymax,betaymin,sbetaymin}=FindMaxMin[betay];

headingbeta = {{"Max[\!\( \[Beta]\_x \)], Max[\!\( \[Beta]\_y \)] [m]",
                " s-locations of Maxs [mm]",
                "Min[\!\( \[Beta]\_x \)], Min[\!\( \[Beta]\_y \)] [m]",
                " s-locations of Mins [mm]"},
               {"x-Horizontal","y-Vertical"}};

outbeta = {
{NumberForm[betaxmax,ndigits],NumberForm[betaymax,ndigits]},
{NumberForm[1000*sbetaxmax,ndigits],NumberForm[1000*sbetaymax,ndigits]},
{NumberForm[betaxmin,ndigits],NumberForm[betaymin,ndigits]},
{NumberForm[1000*sbetaxmin,ndigits],NumberForm[1000*sbetaymin,ndigits]}};

StylePrint[ TableForm[outbeta,TableHeadings -> headingbeta] ];,
(*else*)
Null
];

(* --- Courant-Snyder functions alpha and gamma *)

alphax[s_] := -betax'[s]/2;
alphay[s_] := -betay'[s]/2;

gammax[s_] := (1 + (alphax[s])^2 )/betax[s];
gammay[s_] := (1 + (alphay[s])^2 )/betay[s];

If[OutputAlpha,
  
PlotAlpha = 
Plot[ {alphax[s],alphay[s]}, {s,0,lperiod},
  PlotRange        -> All,
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> {"s (m)",
   "\!\( \[Alpha]\_x \), \!\( \[Alpha]\_y \) \!\( [1] \)",
                 "Undepressed \[Alpha] Functions (black = x, red = y)",None}, 
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle  
    ];

(* --- output alpha plot in eps format *)
Export["im_plt_lat_alpha.eps",PlotAlpha,"EPS"];

(* --- show alpha plot in notebook *)
Print[PlotAlpha] 
];

If[OutputGamma,
  
PlotGamma = 
Plot[ {gammax[s],gammay[s]}, {s,0,lperiod},
  PlotRange        -> All,
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> {"s (m)",
   "\!\( \[Gamma]\_x \), \!\( \[Gamma]\_y \) \!\( [1/Sqrt[m]] \)",
                 "Undepressed \[Gamma] Functions (black = x, red = y)",None}, 
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle  
    ];

(* --- output gamma plot in eps format *)
Export["im_plt_lat_gamma.eps",PlotGamma,"EPS"];
    
(* --- show alpha plot in notebook *)
Print[PlotGamma]
];  
  
(* --- w function *)
    
wx[s_] := Sqrt[betax[s]];
wy[s_] := Sqrt[betay[s]];

If[OutputW, 

PlotW = 
Plot[ {wx[s],wy[s]}, {s,0,lperiod},
  PlotRange        -> All, 
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> {"s (m)",
   "\!\( w\_x \), \!\( w\_y \) \!\( [Sqrt[m]] \)",
                 "Undepressed w Functions (black = x, red = y)",None}, 
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle  
    ];

(* --- output w plot in eps format *)
Export["im_plt_lat_w.eps",PlotW,"EPS"];
    
(* --- show w plot in notebook *)
Print[PlotW] 
];  

(*
Output properties of beam 
*)

(* --- header *)
StylePrint["Beam Properties","Section"];

(* --- information printouts *)
headingbeam = {
{"Dimensionless Perveance, Q",
 "RMS Emittances [mm-mrad]:", 
 "  \!\( \[CurlyEpsilon]\_\( rms,x\) \)",
 "  \!\( \[CurlyEpsilon]\_\( rms,y\) \)", 
 "RMS Edge Emittances [mm-mrad]:",
 "  \!\( \[CurlyEpsilon]\_x = 4*\[CurlyEpsilon]\_\( rms,x\) \)",
 "  \!\( \[CurlyEpsilon]\_y = 4*\[CurlyEpsilon]\_\( rms,y\) \)",
 "Depressed Phase Advances [deg/period]",
 "   x-plane, \!\( \[Sigma]\_x \) [deg/period]",
 "   y-plane, \!\( \[Sigma]\_y \) [deg/period]",
 "Tune Depressions:",
 "  \!\( \[Sigma]\_x /\[Sigma]\_\( 0x \) \)",
 "  \!\( \[Sigma]\_y /\[Sigma]\_\( 0y \) \)"
},None        };

outputbeam = {
 ScientificForm[NumberForm[N[Q],ndigits]],
 " ", 
 NumberForm[NumberForm[N[10^6*emitx/4],ndigits]],
 NumberForm[NumberForm[N[10^6*emity/4],ndigits]],
 " ", 
 NumberForm[NumberForm[N[10^6*emitx],ndigits]],
 NumberForm[NumberForm[N[10^6*emity],ndigits]],
 " ",
 NumberForm[N[sigmax*180/Pi],ndigits],
 NumberForm[N[sigmay*180/Pi],ndigits],
 " ",
 NumberForm[N[sigmax/sigma0x],ndigits],
 NumberForm[N[sigmay/sigma0y],ndigits]
             };

StylePrint[ TableForm[outputbeam,TableHeadings -> headingbeam] ];


(*
Matched Envelope
*)

(* --- header *)
StylePrint["Matched Solution","Section"];


(* --- make plot of matched envelope *)
envscale = If[OutputEnvForm===Edge,1,(* rms *) 1/2];
PlotEnvelope = 
Plot[ {envscale*rx[s]*1000, envscale*ry[s]*1000}, {s,0,lperiod},
  PlotRange        -> All, 
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> If[OutputEnvForm === Edge, 
   {"s [m]","\!\( r\_x \), \!\( r\_y \) [mm]",
    "Matched Envelope rms Edge Radii (black = x, red = y)",None}, 
                        (* else; rms *) 
   {"s [m]","\!\( x\_\(rms\) \), \!\( y\_\(rms\) \) [mm]",
    "Matched Envelope rms Radii (black = x, red = y)",None}
                        ],
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle 
    ];

(* --- output matched envelope plot in eps format *)
Export["im_plt_env_coord.eps",PlotEnvelope,"EPS"];

(* --- show matched envelope plot in notebook *)
Print[PlotEnvelope]; 


(* --- plot of matched envelope angles *)
PlotEnvelopeAngle = 
Plot[ {envscale*rxp[s]*1000, envscale*ryp[s]*1000}, {s,0,lperiod},
  PlotRange        -> All, 
  Axes             -> False, 
  Frame            -> True,
  PlotRangePadding -> {None,Automatic},
  PlotStyle        -> {Black,Red}, 
  FrameLabel       -> If[OutputEnvForm === Edge, 
    {"s [m]","\!\( r\_x\^' \), \!\( r\_y\^' \) [mrad]",
     "Matched Envelope rms Edge Angles (black = x, red = y)",None}, 
                         (* else; rms *)
    {"s [m]","\!\( x\_\(rms\)\^' \), \!\( y\_\(rms\)\^' \) [mrad]",
     "Matched Envelope rms Angles (black = x, red = y)",None}
                        ],
  ImageSize        -> plotwidth, 
  BaseStyle        -> basestyle  
    ];

(* --- output matched envelope angle plot in eps format *)
Export["im_plt_env_ang.eps",PlotEnvelopeAngle,"EPS"];


(* --- show matched envelope angle plot in notebook *)
Print[PlotEnvelopeAngle];


(* --- calculate envelope averages *) 
rxbar   = Bar[rx[s]];
rybar   = Bar[ry[s]];
rxrybar = Bar[rx[s]*ry[s]];

(* --- Calculate envelope max and mins
       This is done in a crude but robust manner to cover all possibilities.
       First the lattice period is divided into small subintervals.  Then
       a numerical minimization is used over the subinterval using the
       interval with the highest/lowest values to obtain a more accurate
       results.  This method does not rely on symmetries and cannot deal with 
       degeneracies (if they exist).   

       Mathematica's optimization routine returns an error if it stops
       on interval endpoints.  So periodic extensions have been defined 
       to allow more robust searches incase of extrema occuring at endpoints.  

       Finally, checks are made to map results back to the fundamental 
       lattice period.   
*)


If[ Not[LatticeType === Continuous],

  (* --- standard case of periodic lattice *)
  {rxmax,srxmax,rxmin,srxmin} = FindMaxMin[rx];
  {rymax,srymax,rymin,srymin} = FindMaxMin[ry];
  {rxpmax,srxpmax,rxpmin,srxpmin} = FindMaxMin[rxp];
  {rypmax,srypmax,rypmin,srypmin} = FindMaxMin[ryp];,
(* --- else, Continuous Focusing Values *) 

  (* --- set "extremum" values to average values for continuous focusing *)
  rxmax = rxbar;
  rxmin = rxbar;
  rymax = rybar;
  rymin = rybar;
  rxpmax = 0.;
  rxpmin = 0.;
  rypmax = 0.;
  rypmin = 0.

  ];

(* --- output envelope max and mins.  *)
If[ OutputEnvForm === Edge, 
headingmatchenv = {
{"Edge Radii, \!\( r\_x = 2\[LeftAngleBracket] x\^2 \[RightAngleBracket]\^\(1/2\),
  r\_y = 2\[LeftAngleBracket] y\^2 \[RightAngleBracket]\^\(1/2\) \) :",
 "  Avg (Lattice Period), \!\( \( r\_x \)\&_, \( r\_y \)\&_ \) [mm]",
 "  Max, \!\( Max[r\_x], Max[r\_y] \) [mm]",
 "    s-locations of Maxs [mm]",
 "  Min, \!\( Min[r\_x], Min[r\_y] \) [mm]",
 "    s-locations of Mins [mm]",
 "Angles, \!\( r\_x\^\('\), r\_y\^\('\) \) :",
 "  Max, \!\( Max[r\_x\^\('\)], Max[r\_y\^\('\)] \) [mrad]",
 "    s-locations of Maxs [mm]",
 "  Min, \!\( Min[r\_x\^\('\)], Min[r\_y\^\('\)] \) [mrad]",
 "    s-locations of Mins [mm]",
 "Matching Conditions:",
 "  Radii,  \!\( r\_x[0], r\_y[0] \) [mm]",
 "  Angles, \!\( r\_x\ ' [0], r\_y ' [0] \) [mrad]"
},{"x-Horizontal","y-Vertical"}
                  },
(* else, rms format *)
headingmatchenv = {
{"Edge Radii, \!\( x\_\(rms\) = \[LeftAngleBracket] x\^2 \[RightAngleBracket]\^\(1/2\),
  y\_\(rms\) = \[LeftAngleBracket] y\^2 \[RightAngleBracket]\^\(1/2\) \) :",
 "  Avg (Lattice Period), \!\( \( x\_\(rms\) \)\&_, \( y\_\(rms\) \)\&_ \) [mm]",
 "  Max, \!\( Max[x\_\(rms\)], Max[y\_\(rms\)] \) [mm]",
 "    s-locations of Maxs [mm]",
 "  Min, \!\( Min[x\_\(rms\)], Min[y\_\(rms\)] \) [mm]",
 "    s-locations of Mins [mm]",
 "Angles, \!\( x\_\(rms\)\^\('\), y\_\(rms\)\^\('\) \) :",
 "  Max, \!\( Max[x\_\(rms\)\^\('\)], Max[y\_\(rms\)\^\('\)] \) [mrad]",
 "    s-locations of Maxs [mm]",
 "  Min, \!\( Min[x\_\(rms\)\^\('\)], Min[y\_\(rms\)\^\('\)] \) [mrad]",
 "    s-locations of Mins [mm]",
 "Matching Conditions:",
 "  Radii,  \!\( x\_\(rms\)[0], y\_\(rms\)[0] \) [mm]",
 "  Angles, \!\( x\_\(rms\)\ ' [0], y\_\(rms\) ' [0] \) [mrad]"
},{"x-Horizontal","y-Vertical"}
                  }
  ];

If[ OutputEnvForm === Edge, 
outputmatchenv = {
 {" ", " "},
 {NumberForm[rxbar*1000,ndigits],NumberForm[rybar*1000,ndigits]},
 {NumberForm[rxmax*1000,ndigits],NumberForm[rymax*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxmax*1000,ndigits],
  NumberForm[srymax*1000,ndigits]}
   ],
 {NumberForm[rxmin*1000,ndigits],NumberForm[rymin*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxmin*1000,ndigits],
  NumberForm[srymin*1000,ndigits]}
   ],
 {" ", " "},
 {NumberForm[rxpmax*1000,ndigits],NumberForm[rypmax*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxpmax*1000,ndigits],
  NumberForm[srypmax*1000,ndigits]}
   ],
 {NumberForm[rxpmin*1000,ndigits],NumberForm[rypmin*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxpmin*1000,ndigits],
  NumberForm[srypmin*1000,ndigits]}
   ],
 {" ", " "},
 {NumberForm[ rx[0]*1000,ndigits],NumberForm[ ry[0]*1000,ndigits]},
 {NumberForm[rxp[0]*1000,ndigits],NumberForm[ryp[0]*1000,ndigits]}
                 }, 
(* else, rms format *)
outputmatchenv = {
 {" ", " "},
 {NumberForm[rxbar/2*1000,ndigits],NumberForm[rybar/2*1000,ndigits]},
 {NumberForm[rxmax/2*1000,ndigits],NumberForm[rymax/2*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxmax*1000,ndigits],
  NumberForm[srymax*1000,ndigits]}
   ],
 {NumberForm[rxmin/2*1000,ndigits],NumberForm[rymin/2*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxmin*1000,ndigits],
  NumberForm[srymin*1000,ndigits]}
   ],
 {" ", " "},
 {NumberForm[rxpmax/2*1000,ndigits],NumberForm[rypmax/2*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxpmax*1000,ndigits],
  NumberForm[srypmax*1000,ndigits]}
   ],
 {NumberForm[rxpmin/2*1000,ndigits],NumberForm[rypmin/2*1000,ndigits]},
 If[LatticeType===Continuous, {"Arb","Arb"}, 
 {NumberForm[srxpmin*1000,ndigits],
  NumberForm[srypmin*1000,ndigits]}
   ],
 {" ", " "},
 {NumberForm[ rx[0]/2*1000,ndigits],NumberForm[ ry[0]/2*1000,ndigits]},
 {NumberForm[rxp[0]/2*1000,ndigits],NumberForm[ryp[0]/2*1000,ndigits]}
                 }
  ];

StylePrint[ TableForm[outputmatchenv,TableHeadings -> headingmatchenv] ];

If[ OutputEnvForm === Edge, 
headingavgrad = {
{"Average Radius Measures:",
 "  \!\( \@\( \(r\_x r\_y\)\&_ \) \) [mm]",
 "  \!\( (\(r\_x\)\&_ + \(r\_y\)\&_ )/2 \) [mm]"
},None
                }, 
(* else: rms form *)
headingavgrad = {
{"Average Radius Measures:",
 "  \!\( \@\( \(x\_\(rms\) y\_\(rms\)\)\&_ \) \) [mm]",
 "  \!\( (\(x\_\(rms\)\)\&_ + \(y\_\(rms\)\)\&_ )/2 \) [mm]"
},None
                }
  ];

If[ OutputEnvForm === Edge, 
outputavgrad = {
 " ",
 NumberForm[Sqrt[rxrybar]*1000,  ndigits],
 NumberForm[(rxbar+rybar)/2*1000,ndigits]
               }, 
(* else: rms form *)
outputavgrad = {
 " ",
 NumberForm[Sqrt[rxrybar]/2*1000,  ndigits],
 NumberForm[(rxbar/2+rybar/2)/2*1000,ndigits]
               }
  ];

StylePrint[ TableForm[outputavgrad,TableHeadings -> headingavgrad] ];

(*
Numerical Parameters  
*)

(* --- header *)
StylePrint["Matched Solution -- Numerical Parameters","Section"];

(* --- information printouts *)
headingmatchnum = {
{"Parameterization Case",
 If[SolCase === 0,"  Hybrid Case"," "], 
 "Specified Fractional Tolerance", 
 "Achieved  Fractional Tolerance",
 "Iterations Needed", 
 "CPU Time for Solution [sec]"
},None            };

outputmatchnum = {
 SolCase,
 If[SolCase === 0,HybridCase," "],
 ScientificForm[N[tol]],
 ScientificForm[NumberForm[N[tolachieved],ndigits]], 
 iterations, 
 time
                 };

StylePrint[ TableForm[outputmatchnum,TableHeadings -> headingmatchnum] ];
