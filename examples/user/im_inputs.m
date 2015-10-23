(* 
im_inputs.m    

User inputs for Iterated Matching (IM) Method program for generating 
matched beam envelopes to to KV envelope equations.  

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

(* ***************************************************************** *)
(* Input Parameters                                                  *)
(* ***************************************************************** *)
(*
The parameters below are employed to define the focusing lattice 
and the beam/matched envelope.  A control variable SolCase is used to 
specify which parameters are employed in construction of the matched 
envelope solution.  If a parameter is not needed for the SolCase 
specified, then the parameter is ignored in constructing the matched 
solution.  Then after construction of the matched solution that value 
of that irrelevant parameter will be reset to a value consistent with 
the matched envelope (specified by the relevant parameters). 

Comments: 

1/ 
We have set up the lattice to work for 
sigma0x = sigma0y.  However, we use two input variables 
in this case as well as retaining the more general form in the 
programs.  This enables more straightforward modification of the code 
to problems with lesser degrees of symmetry.  

2/ 
If high tolerance solutions are constructed, then parameters should be put in
with high precision (infinite, if possible) to ensure that they will not 
limit the accuracy of various calculations carried out in setting up the 
numerical problem.  For example, use 
   
      sigmax = (1/4)*sigma0x;   NOT   sigmax = 0.25*sigma0x;
            -or- 
      sigmax = Rationalize[0.25*sigma0x,10^{-4*$MachinePrecision}]
 
3/ 
In paralleling the paper's notation, at the end of variable names we will 
write all subscripts first and superscripts second.  

  wp .... working precision used in Mathematica

  Solution Case:
  SolCase     .... 0  lattice params + Q, emitx, emity   specified 
                   1  lattice params + Q, TuneDepx, TuneDepy specified 
                   2  lattice params + emitx = emity, TuneDepx = TuneDepy 
                          specified
                      
                      comments: 
                        Case 0 can take much longer to run than Case 1 and 2 
                          since the IM method is more efficient with specified 
                          phase-advances. So this requires much more numerical 
                          work with a relatively slow (Mathematica-based) tools.
                         
                        Case 2 is restricted from the most general 
                          situation to keep the programs in simple form.

                       

  HybridCase  .... 1  For SolCase = 0, use hybrid case 1 with root finding 
                   2  For SolCase = 0, use hybrid case 2 with root finding

  Lattice:

  LatticeType .... flag for lattice type: 
                     Continuous, Solenoid, or Quadrupole for piecewise 
                     constant cases covered in paper, or 
                     UserInput for user specified
  lperiod     .... lattice period [meters] defined on s = [0,lperiod] 
                     Any user input lattice kappax and kappay 
                     will be repeated over this period

  For LatticeType = Solenoid or Quadrupole (ignored otherwise):
    eta         .... fractional focusing occupancy in period,  0 < eta <= 1 
  For LatticeType = Quadrupole 
    alpha       .... lattice syncopation parameter, 0 <= alpha <= 1 

  kappax      .... x-plane linear focusing strength [meters^-2]
  kappay      .... y-plane linear focusing strength [meters^-2]

  kappaxinput[s] .... For LatticeType = UserInput; x-plane linear focusing 
                      function defined s in [0,lperiod]. 
  kappayinput[s] .... For LatticeType = UserInput; y-plane linear focusing 
                      function defined s in [0,lperiod].  

  brkptsinputx,y .... List of "breakpoints" for LatticeType = UserInput
                      in s = [0,lperiod].  The list should include any 
                      discontinuities in piecewise constant lattices.  
                      If none, set as {0,lperiod}.   

  sigma0x     .... x-plane undepressed phase advance [rad/period]
  sigma0y     .... y-plane undepressed phase advance [rad/period]
                     Note: sigma0x and sigma0y *must* be set consistent 
                           with kappax and kappay.  For LatticeType not 
                           UserInput this is done automatically.  For 
                           LatticeType = UserInput sigma0x and sigma0y are 
                           set from kappax and kappay so kappax and kappay 
                           must be set consistent with desired lattice 
                           focusing strength.  
  Beam:

  Q       .... dimensionless perveance 
  emitx   .... x-plane beam rms edge emittance [meter-radian]
  emity   .... y-plane beam rms edge emittance [meter-radian]

  sigmax  .... x-plane depressed phase advance [rad/period] (auto-set) 
  sigmay  .... y-plane depressed phase advance [rad/period] (auto-set) 

  TuneDepx = sigmax/sigma0x  (used input <=> sigmax)
  TuneDepy = sigmay/sigma0y  (used input <=> sigmay) 

  Note: Due to program structure issues TuneDepx, TuneDepy are used 
  in place of sigmax and sigmay as input parameters.  sigmax is calculated 
  from TuneDepx as sigmax = TuneDepx*sigma0x etc.  

*)
(* --- store machine working precision for input of decimal numbers *)
wp   = $MachinePrecision;

(* --- Solution Case *) 
SolCase = 1;
HybridCase = 1;

(* --- Lattice Parameters: 
         use: quantity = Rationalize[#.##,10^(-4*wp)] for decimal numbers *)

LatticeType = UserInput;
lperiod = 1/2;

eta     = 1/2;
alpha   = 1/2;

sigma0x = 80*(Pi/180);
sigma0y = sigma0x;

(* --- --- UserInput Lattice specification defined on [0,lperiod].
           Used only when LatticeType = UserInput and otherwise ignored.  

           Comments: 
             Replace this function with the desired (numerical: must return 
               number) function defined on [0,lperiod] 

             eta, alpha, and sigma0x, sigma0y are ignored when this function 
               is used with LatticeType = UserInput. 

             Example provided below for a sinusoidally varying quadrupole 
               type lattice has amplitude set to achieve 
               simga0x = simga0y = sigma0 = 80 deg/period.  Phase advances 
               are calculated consitent with input amplitudes.   
*)

kappaxinput[s_] := 45.8528*Sin[2*Pi*s/lperiod];
kappayinput[s_] := -kappaxinput[s];

(* --- --- List of break points of any lattice discontinuities in 
           UserInput lattice functions kappaxinput[s] and kappayinput[s] 
           defined on [0,lperiod]. If there are no discontinuities 
           (continuous function), then set to {0,lperiod}.  
*)
 
brkptsinputx = {0,lperiod};
brkptsinputy = {0,lperiod};


(* --- Beam Parameters:
         use: quantity = Rationalize[#.##,10^(-4*wp)] for decimal numbers *)

Q     = 1.*10^(-4);

emitx = 10*10^(-6);
emity = emitx;


TuneDepx = 0.2;
TuneDepy = TuneDepx;


(* ***************************************************************** *)
(* Optional Program Outputs:                                         *)
(* Focus Lattice, Characteristic Orbits, and Envelope Modes          *)
(* ***************************************************************** *)
(* 
Control flags for optional program outputs 
   
 Focuing Lattice 
   OutputBeta   = True   =>  Plot the undepressed betatron functions 
                             describing the applied focusing lattice.  

   OutputAlpha  = True   =>  Plot undepressed alpha functions 
                             describing the applied focusing lattice.  

   OutputGamma  = True   =>  Plot undepressed gamma functions 
                             describing the applied focusing lattice.  
 
   OutputW      = True   =>  Plot undepressed w functions
                             describing the applied focusing lattice. 

 Characteristic Particle Orbits 
   OrbitSolve   = True   =>  Calculate and output characteristic 
                             depressed and undepressed principal orbit 
                             within matched beam envelope.
   
 Envelope Stability Properties 
   EnvStability = True   =>  Calculate and output envelope mode 
                             stability properties using the formulation 
                             by Lund and Bukh, "Stability properties of the 
                             transverse envelope equations describing intense 
                             ion beam transport," PRSTAB 7, 024801 (2004)
*)

OutputBeta   = True;
  
OutputAlpha  = False;
OutputGamma  = False;
OutputW      = False;
  
OrbitSolve   = True;
EnvStability = True;

(* Control parameters used to set both undepressed and depressed oscillations

  nperiod  = number of periods to advance, need not be interger 

  soi      = initial s of orbit, m 
  sof      = final   s of orbit, m   * set from soi and nperiod * 


  fenvic   = fraction edge of initial radial coordinate 
  fedge    = fraction beam edge of particle oscillation 
*)


nperiod = 2*Pi/sigmax;    (* Note: one betatron osc ... may not be integer *)
(* nperiod = 2*Pi/sigma0x; *)

soi = 0;                               
sof = Period;            

xoi  = Fenvx;
xpoi = Femitx; 

yoi  = Fenvy;
ypoi = Femity; 

fenvx   = 0.6;    
femitx  = 0.6;    

fenvy   = 0.6;    
femity  = 0.6;    

(* ***************************************************************** *)
(* Numerical Control Parameters                                      *)
(* ***************************************************************** *)
(*
Parameters for numerical and output control. 
*)


(* --- solution 
   tol     .... fractional accuracy tolerance of matched solution 
   ssample .... list of s values in the lattice period s = [0,lperiod] 
                  to calculate error tolerance 
*)
tol     = 10^(-6);
ssample = lperiod*Table[i,{i,0,100}]/100;


(* --- Set up initial conditions for principal orbits. 
       These apply to both the seed iteration and the general iteration. 

   si .... initial s coordinate s_i to apply initial conditions 

   cxi  .... c_x(s_i)
   cxpi .... d c_x (s)/ds |_{s=s_i}

   etc.  
*)
si   = 0;
 
cxi  = 1;  cxpi = 0;
sxi  = 0;  sxpi = 1;

cyi  = 1;  cypi = 0; 
syi  = 0;  sypi = 1;


(* ---- Numerical control parameters.  

   Accuracies and precisions are specified as x => 10^(-x)

   For numerical integration of principal particle orbits:
   ndag     .... numerical differential equations, accuracy goal 
   ndpg     .... numerical differential equations, precision goal 
   ndmethod .... numerical method specification:
                   Automatic, ExplicitRungeKutta, 
                   {"SymplecticPartitionedRungeKutta",
                    "DifferenceOrder"   -> Automatic,
                    "PositionVariables" -> {x[s]} 
                   } 
     Method comments:
       - Orbits are only integrated for one lattice period which, for low 
         accuracy results, makes the choice of method less critical.  
       - Automatic: Default method appear to be fast and switches between 
         several methods and orders as dictated by accuracy needs. 
       - ExplicitRungaKutta: Less efficient but removes ambiguity in method 
         properties.  Variable order.  
       - SymplecticPartitionedRungeKutta: RungaKutta method modified to 
         preserve Hamiltonian structure of dynamics.  Slower, but a good choice
         for high accuracy applications. 
  
   For constraint integrals of trial envelope functions:
   niag     .... numerical integrals, accuracy  goal 
   nipg     .... numerical integrals, precision goal 
   nimethod .... numerical method specification:
                   Automatic  

   For root finding to implement Case 0 hybrid solutions:
   rfag     .... root finding, accuracy  goal 
   rfpg     .... root finding, precision goal 
   rfmethod .... numerical method specification:
                   Automatic, Newton, Secant
    Method comments:
      - The settings of tol and the root finding accuracy (rfag and rfpg) WILL 
        interact.  For very high accuracy root finds, care should be taken to 
        set tol for higher accuracy than specified in the root finding to 
        avoid noise issues.  Numerical accuracies of differential equation 
        solves and integrations should also be set consistently.   
      - Two dimensional root finding is carried out with only 
        function evaluations rather than function and Jacobian evaluations.  
        Derivatives must be found by finite differencing.  Logically, in such 
        cases, 2-point starting initial values for the root finding 
        should be given.   
      - Automatic: Default, appears to be a Newton-like method with good 
        performance. 
      - Newton: In Mathematica 5.2 returns a warning with 2-point initial 
        starting values saying that the 2nd point is ignored.  Since the 
        program cannot calculate analytical Jacobian values from the function 
        the root finding is applied to, this means that the root finding 
        program must take extra steps to determine the scale to first finite 
        difference the function.  Hopefully this will be corrected in later 
        versions of Mathematica.   The program appears to determine the needed 
        scale efficiently rendering this (inappropriate) program structure of 
        little practical consequence.  
      - Secant: Appears to employ a Secant like multi-dimensional 
        generalization of the standard 1d Secant method.  The Secant 
        method may preform better than the Newton method when 
        there is numerical noise in the root finding functions.  


   General Comments:
     - Original values of numerical parameters distributed in this file are 
       set for reasonably fast running and to be consistent with the specified 
       fractional accuracy tol.   
     - If higher accuracy (smaller tol values) are desired, numerical 
       parameters and methods must be set consistently and solutions tested 
       to verify accuracies.   
     - Generally set AccuracyGoal variables (xxag) to several digits less 
       than the machine working precision to deemphasize absolute accuracy 
       measures and use relative precisions.  Setting to Infinity can cause 
       problems though.  
     - Generally set PrecisionGoal variables (xxpg) to values corresponding 
       to MORE precision than the specified tol.  Some guidance:  
         ndpg  .... particle orbits, set most accurate 
         nipg  .... integrals of envelopes, set a little less accurate than 
                    orbits. 
         rfpg  .... used only for SolCase = 0 for root finding.  Set least 
                    accurate. 
  Comments:
    ndag, ndpg:  Had problems with NDSolve when ndag was changed from default 
                 value Automatic when ndpg was also set.  Setting 
                 ndag = Automatic, and ndpg = Infinity (to not consider) 
                 works.  
*) 
ndag = Automatic;
ndpg = Infinity;
ndmethod = Automatic;

niag = wp - 2;
nipg = 8;
nimethod = Automatic;

rfag = wp - 4;
rfpg = 4;
rfmethod = Automatic;

(* --- debugging printouts and plots 
         See Also: DeBug flag set at top of im_solver.m for further debug info.  
*)
IterationPrint = False;      (* Set True or False *)
IterationPlot  = False;      (* Set True or False *)

(* --- output control *)
ndigits = 5;                        (* Number significant figures for various 
                                         output quantities *)
basestyle   = {"Times-Roman",14};   (* Font for plot labels *)
plotwidth   = 500;                  (* plot width in notebook  *)


