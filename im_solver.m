(* 
im_solver.m    

Main control program for Iterated Matching (IM) Method for constructing 
matched beam envelopes to to KV envelope equations.  

See readme.txt for an overview of program operation.  The program can be 
executed by loading im_solver.m within a notebook:
   >> im_solver.m [shift-return]
Other than to change the program, only system and numerical parameters and 
program control flags contained with im_inputs.m need be changed.   

Steve Lund and Sven Chilton 
Lawrence Berkeley National Lab
Jan, 2006

Updated by Steve Lund and Kei Fukushima (Hiroshima University) 
June, 2012  
Updated by Steve Lund 
March, 2014; December 2016   

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
Debug Printout Flag: Use to help track down code issues.  
  See Also: IterationPrint and IterationPlot in im_inputs.m  
*)
DeBug = True;
DeBug = False;
  
(* ***************************************************************** *)
(* Version Check                                                     *)
(* ***************************************************************** *)
(* 
Snytax changes in Mathematica require version 5.0 or higher for 
correct program operation.  If desired, it should be possible to 
easily modify the program run in earlier versions by changing 
a few FindRoot[] call arguments to reflect older syntax.  
*)

If[ DeBug || $VersionNumber < 8.0,
    Print["Program written for Mathematica 8.0 or higher."];
    Print["Verified to work with Mathematica 10.0 3/2015."];
    Print["May work for versions 5.0 and older with minor modifications."];
    If[ $VersionNumber < 8.0, Abort[]]
  ];


(* ***************************************************************** *)
(* Input Parameters                                                  *)
(* ***************************************************************** *)
(* 
Read in system and numerical parameters plus program output control flags.  
This should contain the full user description of the lattice and 
system parameters.  
*)

If[DeBug, Print["DeBug: Import im_inputs.m"]];
<< im_inputs.m 
(* Abort[]; *)   

(* **************************************************************** *)
(* Name Notebook and Start Timing                                   *)
(* **************************************************************** *)

(* --- Name the notebook input.  The save command can only be read in once.
 
       Comment:  
       Was seeing odd bug. My install of mathematica is trying 
       to find a deleted directory to save the notebook to.  I have no idea 
       why this is happening. Perhaps there is some directory info saved in 
       the program preferences.  I don't think this should impact others.  
       But it will result in the notebook not being named and saved.  

*)

If[$Notebooks,
  MatchEnvSolve = SelectedNotebook[ ];
  If[Not[ValueQ[nbsave]],NotebookSave[MatchEnvSolve,"MatchEnvSolve"];
                         nbsave=False
    ]
  ];

(* --- start timing to measure cpu time in seconds *)
timestart = TimeUsed[];

 
(* **************************************************************** *)
(* Output Header                                                    *)
(* **************************************************************** *)
(*
Generate and print header with date, user, and machine information. 
*)
 
Print[Style["Matched Envelope Solution -- IM Method","Subtitle"]];
rundate = Date[];
Print[Style[ToString[Part[rundate,2]]<>"-"<>ToString[Part[rundate,3]]<>
           "-"<>ToString[Part[rundate,1]]<>"   by "<>ToString[$UserName]<>
           " on "<>ToString[$MachineName],"Subsubtitle"]
     ];
Print[Style["Code Provided by Steve Lund","Subsubtitle"]];
Print[Style["Michigan State University (MSU), Facility for Rare Isotope Beams (FRIB)","Subsubtitle"]];


(* ***************************************************************** *)
(* Setup Utilities                                                   *)
(* ***************************************************************** *)
(*
Read in basic utility functions employed in the IM method 
*)
If[DeBug, Print["DeBug: Import im_utilities.m"]];  
<< "im_utilities.m"
(* Abort[]; *)  
  
(* ***************************************************************** *)
(* Applied Focusing Lattice                                          *)
(* ***************************************************************** *)
(*
Generate linear applied focusing lattice 
*)

<< "im_lattice.m" 
If[DeBug, Print["DeBug: Read in im_lattice.m"]];
(* Abort[]; *)
  
(* ***************************************************************** *)
(* Partial Error Checks                                              *)
(* ***************************************************************** *)
(*
Read in a partial set of error traps  
*)
If[DeBug, Print["DeBug: Import im_err_chk.m"]];
<< "im_err_chk.m" 
(* Abort[]; *)
  
(* ***************************************************************** *)
(* Calculate Matched Envelope to Tolerance                           *)
(* ***************************************************************** *)

(* Set flag to run IMSolver for standard matched envelope cases *)
IMSolve = True;

(* --- Direct solution case for any initial condition.
         This can produce a potentially mismatched envelope 
*) 

If[DeBug, Print["DeBug: Generating direct solution if applicable."]];
If[ SolCase == - 1,
  solution = DirectSolve[];
  (* *)
  rx  = xrx  /. solution[[1]][[1]];
  rxp = xrxp /. solution[[1]][[2]];
  ry  = xry  /. solution[[1]][[3]];
  ryp = xryp /. solution[[1]][[4]];
  (* *)
  tolachieved = "NA";
  iterations  = "NA"; 
  IMSolve = False 
  ];
(* Abort[]; *)

(* --- read in functions to construct continuous limit envelopes used 
        in forming a seed iteration 
*)
If[DeBug, Print["DeBug: Import im_cont.m"]];
If[IMSolve, << "im_cont.m" ];
(* << "im_cont.m" *)
(* Abort[]; *)

(* --- read in functions to generate seed iteration *)
If[DeBug, Print["DeBug: Import im_seed.m"]];
If[IMSolve, << "im_seed.m" ];
(* << "im_seed.m" *)
(* Abort[]; *)

  
(* --- read in functions for general IM method iterations and generation 
       of matched envelope solutions to tolerance  *)
If[DeBug, Print["DeBug: Import im_iterate.m"]];
If[IMSolve, << "im_iterate.m" ];
(* << "im_iterate.m" *) 
(* Abort[]; *)
  
(* --- construct matched envelope solution using IM method; save
         iterations  = iteration count 
         tolachieved = fractional numerical tolerance of solution             
*)
If[DeBug,
Print["DeBug: Calling Match[] to generate matched envelope"];
Print["       Set IterationPrint = True in im_inputs.m for mode info."];
Print["       Set IterationPlot  = True in im_inputs.m for mode info."]
  ];     

If[IMSolve, {iterations, tolachieved} = Match[] ];
(* Abort[]; *)

(* --- calculate envelope angles from matched envelope functions  
         only needed in matched iteration case: already calculated in 
         direct integration when SolCase = -1 <=> IMSolve = False
   rxp[s] .... r_x'(s) defined on s = [0,lperiod]
   ryp[s] .... r_y'(s) defined on s = [0,lperiod]
*)
If[ IMSolve, 
  rxp[s_] := rx'[s];
  ryp[s_] := ry'[s]
  ];

(* --- calculate/reset initial conditions 
         Used for SolCase = -1, but reset for matched cases in case user applies.
         Reset will not matter for SolCase = -1 so keep simple. 
*) 
rxi = rx[0];
ryi = ry[0];
rxpi = rxp[0];
rypi = ryp[0];

rxf = rx[lperiod];
ryf = ry[lperiod];
rxpf = rxp[lperiod];
rypf = ryp[lperiod];



(* --- stop timing, so processor time in seconds to generate 
       solution is counted *)
timestop = TimeUsed[];
time     = timestop - timestart;


(* ***************************************************************** *)
(* Diagnostic Output                                                 *)
(* ***************************************************************** *)
(* 
Output formatted information on envelope solution including 
numerical solution of envelope maximums/minimums etc.   

Comment: In case of program issues, the DeBug option should be turned on to 
help locate problems since diagnostic output is ordinarily only made *after* 
the solution is numerically generated.  
*)

If[DeBug, Print["DeBug: Import im_diag.m"]];
<< "im_diag.m" 
(* Abort[]; *)   


(* ***************************************************************** *)
(* Characteristic Particle Orbits                                    *)
(* ***************************************************************** *)
(*
Optional: Make plots of characteristic depressed and undepressed 
particle orbits in matched beam envelope.  
  *)
  
If[DeBug && OrbitSolve, Print["DeBug: Import im_orbit.m"]];
If[ OrbitSolve, << "im_orbit.m"];
(* Abort[]; *) 

(* ***************************************************************** *)
(* Envelope Stability                                                *)
(* ***************************************************************** *)
(*
Optional: Carry out envelope stability analysis of matched beam envelope.  
*)

If[DeBug && EnvStability, Print["DeBug: Import im_env_stability.m"]]; 
If[ EnvStability, << "im_env_stability.m"];
(* Abort[]; *) 
