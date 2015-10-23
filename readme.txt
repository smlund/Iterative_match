**************************** 
Overview   
test mod
****************************
This package of Mathematica programs numerically constructs matched
solutions to the KV envelope equations describing a charged particle beam
in a linear, periodic transport lattice using the "Iterative Matching" 
(IM) method.  The IM formulation is described in:

  "Efficient computation of matched solutions of the KV envelope equations
   for periodic focusing lattices" 

  S.M. Lund^1, S.H. Chilton^2, and E.P. Lee^2

  1. Lawrence Livermore National Laboratory, Livermore, CA 94550 
  2. Lawrence Berkeley  National Laboratory, Berkeley,  CA 94720   

  Physical Review Special Topics -- Accelerators and Beams 
  Vol. 9, 064201 (2006) 

  Paper can be downloaded for free at:
  http://prst-ab.aps.org/

  http://prst-ab.aps.org/abstract/PRSTAB/v9/i6/e064201

The programs both demonstrate the method described in the paper and 
also provide a useful and efficient program for generation of matched 
envelope solutions and analyzing the properties of the solutions.  
Variable names in the programs reflect the notation used in the 
paper.  The paper serves as documentation and should be 
referred to for a detailed understanding of the program.  The program works 
with various common parameterizations of the system that are employed in 
accelerator physics including the common perveance + emittance parameterization
or parameterizations where either the perveance or emittances are replaced 
by phase advance parameters.  The primary limitation of the program is that 
it does not work in the limit of full space-charge depression.  This 
limitation could be removed through simple code extensions, but the authors
decided against implementing this to keep the program relatively clean as is
and this limit is not possible in the laboratory anyway.  
Linear focusing functions corresponding to the piecewise constant lattices 
using in the paper (continuous, solenoidal, and quadrupole doublet) are 
included as options.  Options also allow the user to input any 
x- and y-plane lattice focusing functions desired by simply defining the
desired lattice focusing functions over one lattice period and setting a flag
to tell the program to use them.  An example included 
for a simple quadrupole lattice with sinusoidally varying focusing 
functions is defined in the input file as a representative example which can
be easily modified.  

Please cite the final, published paper if this work/program is 
referenced.  Users are welcome to modify this program for their own 
purposes.  But the authors request that some citation be made on the source 
of material. 

For further information, please contact the primary author:

Steven M. Lund
Facility for Rare Isotope Beams 
Michigan State University 
640 South Shaw Lane 
East Lansing, MI 48824 
 
lund@frib.msu.edu 
517-908-7291 

Comments on improvements are welcome! 

**************************** 
Obtaining the Code   
****************************

The authors originally posted the code arXiv.org linked with a draft of the
PRSTAB manuscript in 2006. This code has been updated to work with version
changes in Mathematica and numerous improvements have been made. The best
place to obtain the code is on the US Particle Accelerator School (USPAS)
course web site for "Beam Physics with Intense Space-Charge" taught
by Lund and Barnard. See:

  https://people.nscl.msu.edu/~lund/uspas/bpisc_2015/

In particular, under lund lectures on "Transverse Centroid and Envelope
Descriptions of Beam Evolution":

  https://people.nscl.msu.edu/~lund/uspas/bpisc_2015/lec_set_07/
      tce.pdf        ... lecture notes
      tce_ho.pdf     ... lecture notes handout formatted (4 slides/page)
      env_match_code ... directory containing program files 

The authors intend to post updates and/or extensions of the code as prove
necessary. If there are problems, you are encouraged to first check for
the most up to date version and then contact us (see info above) if there
are problems.  

The code can also be obtained on github.com.   

Efforts have been made to employ basic programming structures that should 
remain backward compatible when Mathematica is updated. Unfortunately,
Mathematica has had many changes to fairly basic syntax and default
modes of operation which have impacted to the code in subtle ways
as the software has changed version-to-version. In an effort to help
track down issues that may develop with future version changes, the authors
have included debugging provisions. To use these, set in:

  im_solver.m    =>   DeBug = True;   near top of file 
  im_inputs.m    =>   IterationPrint = True; 
                      IterationPlot  = False;

and diagnostic outputs will be made that may assist in identifying problems.   

**************************** 
Running the Code   
****************************
The program can be executed from a Mathematica notebook in 
Linux/Unix by placing all the needed im_*.m files in the directory 
that the notebook is run from, opening a Mathematica notebook, 
and then executing the command:

   << im_solver.m [shift-return]

from within the notebook. The problem setup is defined in:

   im_inputs.m

which is a simple ascii text file containing parameters and flags (Mathematica
syntax) that can be edited with any standard text editors. On some
MAC OSX installations or MS windows, it may be necessary to first execute
the command 

SetDirectory["directory"] 

before executing im_solver.m in the notebook.  Here, "directory" is the 
full path name of the directory containing the *.m program files.    

The file im_inputs.m contains input parameters needed to 
describe the beam, focusing lattice, envelope solution method, and 
numerical and output control parameters. Numerical parameters may need 
to be reset for high precision and special applications. Comments
in im_inputs.m provide guidance on how the numerical parameters 
can be set. Other than to change and update the program, no modifications 
other than to im_inputs.m should be necessary. All the im_*.m files can 
be modified with a simple text editor. Notebooks 
need not be saved since the programs are compact and can be easily 
rerun.  The programs are written to work with the common cases of 
parameterizations of the matched envelope solution as explained in the 
paper (Section II). This results in more elaborate case structures 
being required in the programs.  Various parts of the program such 
as the principal orbit calculations are written for general initial conditions 
to enable simple generalization for other applications.  If desired, a 
shorter version of the program can be constructed for a single 
parameterization case of interest with a specific coordinate system etc, by 
straightforward simplifications of the functions defined.  

**************************** 
Program Files    
****************************

im_solver.tar  .... tar archive file of the full distribution of files listed
                      below.  

readme.txt     .... This file

examples       .... Directory of example runs with output.  Contains
                      sub-directories:
		         sol   ... hard-edge solenoid   lattice example 
			 quad  ... hard-edge quadrupole lattice example
			 user  ... user input lattice example 

im_solver.m    .... Main control program.  This script reads in problem and 
                      description and numerical control parameters contained 
                      in im_inputs.m, then reads in various function 
                      definition files that follow, and then 
                      makes calls to appropriate functions to generate the 
                      matched envelope solution with the IM method.  After
		      generating the solution the 
                      script im_diag.m is run to generate output plots and 
                      information on the matched solution.  If set, optional 
                      particle orbit and envelope stability programs are 
                      executed with more diagnostics.  

im_inputs.m    .... Input parameters for beam description, linear focusing 
                      lattice, and numerical and output control.  All user 
                      inputs are contained in this file.  Arbitrary linear 
                      focusing lattices can be described using Mathematica 
                      functions.  Depending on the run mode, some inputs are 
                      not used and will be reset consistently with the 
                      parameters employed on generation of the matched 
                      envelope.  

im_lattice.m   .... Defines periodic applied focusing lattices, calculates
                       undepressed betatron functions, and calculates 
                       undepressed single particle phase advance inputs
		       (if necessary) for consistency with the lattice.  

im_utilities.m .... Defines various utilities used in parts of the program.  

im_err_chk.m   .... Partial error checks for program inputs.   

im_cont.m      .... Function definitions for the continuous limit which 
                      is used to generate a seed iteration.   

im_seed.m      .... Function definitions to generate a seed (initial) iteration.

im_iterate.m   .... Function definitions to carry out method iterations and 
                      matched envelope solutions under various parameterization
                      cases. This is the main core of the algorithm.   
                   
im_diag.m      .... Generates formatted output information and plots on 
                      the lattice and the matched envelope solution.  Average, 
                      max, and min values of the envelope radii and angles 
                      are calculated.  Output is displayed in the Mathematica 
                      notebook and can be printed from the notebook.  Plots 
                      generated (if turned on in some cases) are also
		      exported as eps files:
		      
                        im_plt_lat_kappa.eps ... lattice focusing functions
			
			im_plt_lat_beta.eps  ... lattice undep betatron funcs
                        im_plt_lat_alpha.eps ... lattice undep alpha funcs
			im_plt_lat_gamma.eps ... lattice undep gamma funcs

                            (beta,alpha,gamma) are the usual Twiss functions
			    specifiying the orientation of the
			    Courant-Snyder invariant ellipse. With beta being
			    the usual undepressed betatron function.   

                        im_plt_env_coord.eps ... matched envelope
                        im_plt_env_ang.eps   ... matched envelope angles


im_orbit.m 
               .... Optional: Calculates depressed and undepressed principal 
                      orbits within the matched beam envelope and outputs 
                      a formated plot. Initial conditions can be set in
		      im_inputs.m and provisions are given to choose initial
		      conditions that remain in the envelope and choose a
		      full oscillation period.  If run, this outputs the
		      plots:

                        im_plt_orbit_x.eps  ... x-plane particle orbit 
			im_plt_orbit_y.eps  ... y-plane particle orbit

                        im_plt_orbit_xp.eps ... x'particle orbit
			im_plt_orbit_yp.eps ... y'particle orbit 

                      Note: if solenoidal focusing is used, these orbits are
		      implicitly in the rotating Larmor frame.  

im_env_stability.m  
               .... Optional: Calculates stability properties of envelope 
                      modes on the matched envelope and outputs formatted 
                      information on the results and a plot of the envelope 
                      eigenvalues to help interpret stability properties.
		      This relies heavily on treatments presented in:
		       Lund and Bukh, "Stability properties of the
		       transverse envelope equations describing intense
		       ion beam transport," Phys. Rev. Special Topics -- 
                       Accelerators and Beams 7, 024801 (2004).  If run, this
		       outputs the plots:

                         im_plt_env_stab_evals.eps  ... envelope eigenvalues

**************************** 
Examples in the distribution  
****************************

Example versions of runs are stored with the source distribution on the USPAS
web site (see above) in a directory "examples" with sub-directories:

  "sol"  solenoid             piecewise continuous focusing lattice 
  "quad" fodo quadrupole      piecewise continuous focusing lattice
  "user" user input focusing with sinusoidally varying kappa functions

In each case, the code (various im_*.m files) used to generate the output, the
Mathematica notebook generated when these are read in, output eps files, and
a pdf export of the notebook are stored. These cases can be rerun to verify
that the code is working correctly.  

All have lattices with 1/2 meter period, sigma0 = 80 degrees undepressed phase
advance per period, are run with SolCase = 1, perveance Q = 1.*10^(-4),
and space-charge tune depression sigma/sigma0 = TuneDep = 0.2
(emittances are adjusted for consistency and output). The focusing
strength (kappa) of the user input quadrupole lattice is input such that
sigma0 = 80 degrees/period. The piecewise constant lattices all have occupancy
eta = 1/2 (period 50% filled with solenoid or quadrupoles respectively). These
examples can be easily modified to run other solution cases (SolCase) for other
parameterizations.   

**************************** 
History  
****************************
The programs were originally written and debugged under Mathematica 
5.1 and 5.2 running on linux/unix platforms during the late 2005 and 
early 2006 calendar years.  A significant update to the program was made 
in 2012 using Mathematica 8.0.  Syntax changes in function argument passing
were made in 2015 so the program would work with Mathematica 10.0 and other
improvements were made. These modifications should be backward compatible.

A full version of this program is implemented in the Warp code for intense
beam simulation. Information on Warp can be found on:

   http://warp.lbl.gov/

Information on this implementation is covered in the Mater's thesis:

  "Implementation of an iterative matching scheme for the
   Kapchinskij-Vladimirskij equations in the WARP code"
   Sven Chilton
   University of California at Berkeley
   Department of Nuclear Engineering
   Spring 2008

