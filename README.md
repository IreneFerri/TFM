# alpha 3 states


************************************************************************
*****  LANDAU WANG ALGORITHM FOR THE ALPHA-3-STATES MODEL **************
************************************************************************

The program has a bug/problem descibed in 'issues', and is not able 
to visit all the energy bins of the histogram for a number of agents 
N > 14

CONTENTS of the 'alpha3states_wl' folder: --------------------------------

 *  'Code' (subfolder):
      'main.f90'
      'module.f90'
      'input.dat' (output name: nom, size: N (must be even), random seed: iseed,
                   initial modifier: ln_f, alpha value: alpha, number of MC attempts: MCtot,
                   flatness tolerance: tol, modifier threshold: eps, numerical tolerance: eps2))
      'r1279.o', 'ran2.o', secs.o' (random generator files)

 *  'Images' (subfolder):
      'ln_g_N10.png', 'ln_g_N12.png', 'ln_g_N14.png' (Plots of the ln(states density(E)) .vs. E)
    
 *  'Scripts'
      'ln_g.gnu' (For the output visualization) Change the filename and 'N = 10' 
                                                for the desired size


COMPILING INSTRUCTIONS (from local folder 'Code'): ------------------------------------------------

gfortran -g -fcheck=all -Wall -fbacktrace module.f90 main.f90 r1279.o ran2.o secs.o -o main.exe

./main.exe (execution)
    
The output is a .csv file with 2 columns. The first one are the energy bins, 
and the second one correspond to the logarithm of the states density corresponding to 
the given energy.

